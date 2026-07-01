/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynamicCode.H"
#include "OSHA1stream.H"
#include "dlLibraryTable.H"
#include "regIOobject.H"
#include "Pstream.H"
#include "stringOps.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"

#include "parse_wmake.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicCode, 0);
}

int Foam::dynamicCode::allowSystemOperations
(
    Foam::debug::infoSwitch("allowSystemOperations", 0)
);

const Foam::fileName Foam::dynamicCode::codeTemplateDirName
    = "codeTemplates/dynamicCode";
#ifndef WM_OPTIONS
#error "WM_OPTIONS must be defined as macro, in compile command"
#else
const char* const wm_options_string = STR(WM_OPTIONS);
#endif
//const char* const Foam::dynamicCode::libTargetRoot =
//    "LIB = $(PWD)/../platforms/$(WM_OPTIONS)/lib/lib";

const Foam::word Foam::dynamicCode::topDirName
(
    "dynamicCode"
);

const char* const Foam::dynamicCode::libTargetRoot
(
    "LIB = $(PWD)/../platforms/$(WM_OPTIONS)/lib/lib"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicCode::addLineDirective
(
    string& code,
    const label lineNum,
    const fileName& name
)
{
    code = "#line " + Foam::name(lineNum) + " \"" + name + "\"\n" + code;
}


void Foam::dynamicCode::copyAndFilter
(
    ISstream& is,
    OSstream& os,
    const HashTable<string>& mapping
)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << "Failed opening for reading " << is.name()
            << exit(FatalError);
    }

    if (!os.good())
    {
        FatalErrorInFunction
            << "Failed writing " << os.name()
            << exit(FatalError);
    }

    // Copy file while rewriting $VARS and ${VARS}
    string line;
    do
    {
        // Read the next line without continuation
        is.getLine(line, false);

        // Expand according to mapping.
        // Expanding according to env variables might cause too many
        // surprises
        stringOps::inplaceExpandCodeTemplate(line, mapping);
        os.writeQuoted(line, false) << nl;
    }
    while (is.good());
}


bool Foam::dynamicCode::resolveTemplates
(
    const wordList& templateNames,
    DynamicList<fileName>& resolvedFiles,
    DynamicList<fileName>& badFiles
)
{
    bool allOkay = true;
    forAll(templateNames, fileI)
    {
        const fileName& templateName = templateNames[fileI];

        const fileName file
        (
            findConfigFile
            (
                templateName,
                dynamicCode::codeTemplateDirName,
                "system"
            )
        );

        if (file.empty())
        {
            badFiles.append(templateName);
            allOkay = false;
        }
        else
        {
            resolvedFiles.append(file);
        }
    }

    return allOkay;
}


bool Foam::dynamicCode::createMakeFiles() const
{
    // Create Make/files
    if (compileFiles_.empty())
    {
        return false;
    }

    const fileName dstFile(this->codePath()/"CMakeLists.txt");

    // Create dir
    mkDir(dstFile.path());

    OFstream os(dstFile);

    if (!os.good())
    {
        FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
    }
    // os  << nl << dynamicCode::libTargetRoot << codeSha1Name_ << nl;

    os<<"cmake_minimum_required(VERSION 3.28)"<<nl;

    os<<"project("<<codeSha1Name_<<" LANGUAGES CXX)"<<nl;

    os<<"set(target_name "<<codeSha1Name_<<")"<<nl;
    //$(PWD)/../platforms/$(WM_OPTIONS)/lib/lib";
    os<<"set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../platforms/"<<wm_options_string<<"/lib)"<<nl;

    os<<"add_library(${target_name} SHARED "<<nl;
    // Write compile files
    forAll(compileFiles_, fileI)
    {
        os<<"    ";
        os.writeQuoted(compileFiles_[fileI], false) << nl;
    }

    os<<")"<<nl;

    os<<"include(options.cmake)"<<nl;

    return true;
}


bool Foam::dynamicCode::createMakeOptions() const
{
    if (compileFiles_.empty())
    {
        return false;
    }

    const fileName dstFile(this->codePath()/"options.cmake");

    // Create dir
    mkDir(dstFile.path());

    OFstream os(dstFile);

    if (!os.good())
    {
        FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
    }

    auto vars = wmakeParse::get_environment_variables();

    wmakeParse::wmake_parse_option option{};
    option.when_undefined_reference=wmakeParse::undefined_reference_behavior::throw_exception;
    std::string combined_options = this->optionsString_+'\n'+this->libsString_;
    // Remove line continuation
    while (true) {
        std::string_view line_continuation{"\\\n"};
        auto pos_beg = combined_options.find(line_continuation);
        if (pos_beg==std::string::npos) {
            break;
        }
        combined_options.replace(pos_beg,line_continuation.size(),"\n");
    }

    auto direct_options = wmakeParse::parse_wmake_file(combined_options ,vars,option);

    std::vector<std::string> link_lib_names;
    std::vector<std::string> include_dirs;
    for (auto it=direct_options.begin();it not_eq direct_options.end();) {
        auto parsed_lib_name = wmakeParse::parse_link_libs(*it);
        if (not parsed_lib_name.empty()) {
            link_lib_names.emplace_back(parsed_lib_name);
            it=direct_options.erase(it);
            continue;
        }
        auto include_dir = wmakeParse::parse_include_dirs(*it);
        if (not include_dir.empty()) {
            include_dirs.emplace_back(include_dir);
            it=direct_options.erase(it);
            continue;
        }
        ++it;
    }

    os<<"# Original value of makeOptions: \n# ";
    for(char ch:combined_options) {
      os<<ch;
      if(ch=='\n') {
        os<<"# ";
      }
    }
    os<<nl<<nl;

    os<<"find_package(Mikeno CONFIG REQUIRED)"<<nl;
    {
        os<<"target_compile_options(${target_name} PRIVATE"<<nl;
        auto it= vars.find("EXE_INC");
        if(it not_eq vars.end()) {
            os<<"    "<<it->second.c_str()<<nl;
        }
        for (const auto & direct_option_str:direct_options) {
            os<<"    "<<direct_option_str.c_str()<<nl;
        }
        os<<")"<<nl;
    }

    {
        os<<"target_include_directories(${target_name} PRIVATE"<<nl;
        for (const auto & dir:include_dirs) {
            os<<"    "<<dir<<nl;
        }
        os<<")"<<nl;
    }

    {
        os<<"set(link_lib_names "<<nl;
        for (auto &lib_name:link_lib_names) {
            os<<"  "<<lib_name<<nl;
        }
        os<<")"<<nl;
        os<<R"(
# Translate -l<> options into cmake. Example: -lsampling -> Mikeno::sampling
foreach (lib_name ${link_lib_names})
    if(TARGET Mikeno::${lib_name})
        target_link_libraries(${target_name} PRIVATE Mikeno::${lib_name})
    else ()
        target_link_libraries(${target_name} PRIVATE ${lib_name})
    endif ()
endforeach ())"<<nl;
    }

    {
        os<<"target_link_libraries(${target_name} PRIVATE"<<nl;
        os<<"\n"
            "    Mikeno::OpenFOAM_Defines\n"
            "    Mikeno::OpenFOAM\n"
            "    Mikeno::OSspecific\n"
            "    Mikeno::finiteVolume\n\n";

        auto it=vars.find("LIB_LIBS");
        if(it not_eq vars.end()) {
            os<<"    "<<it->second.c_str()<<nl;
        }
        os<<")"<<nl;
    }

    return true;
}


bool Foam::dynamicCode::writeDigest() const
{
    const fileName file(digestFile());
    mkDir(file.path());

    OFstream os(file);
    os  << '_';
    os.writeQuoted(sha1_.str(), false) << nl;

    return os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCode::dynamicCode
(
    const dictionary& contextDict,
    const dictionary& codeDict,
    const word& codeName,
    const word& codeDirName,
    const wordList& codeKeys,
    const wordList& codeDictVars,
    const word& optionsFileName,
    const wordList& compileFiles,
    const wordList& copyFiles
)
:
    codeRoot_
    (
        stringOps::expandEnvVar("$FOAM_CASE")/topDirName
    ),
    libSubDir_(stringOps::expandEnvVar("platforms/$WM_OPTIONS/lib")),
    codeName_(codeName),
    codeKeys_(codeKeys),
    codeDictVars_(codeDictVars),
    optionsFileName_(optionsFileName),
    compileFiles_(compileFiles),
    copyFiles_(copyFiles),
    codeStrings_(codeKeys.size())
{
    if (isAdministrator())
    {
        FatalIOErrorInFunction(contextDict)
            << "This code should not be executed by someone with administrator"
            << " rights due to security reasons." << nl
            << "(it writes a shared library which then gets loaded "
            << "using dlopen)"
            << exit(FatalIOError);
    }

    if (!allowSystemOperations)
    {
        FatalIOErrorInFunction(contextDict)
            << "Loading a shared library using case-supplied code is not"
            << " enabled by default" << nl
            << "because of security issues. If you trust the code you can"
            << " enable this" << nl
            << "facility be adding to the InfoSwitches setting in the system"
            << " configDict:" << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system configDict is either" << nl << nl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/configDict" << nl << nl
            << "or" << nl << nl
            << "    $WM_PROJECT_DIR/etc/configDict" << nl
            << endl
            << exit(FatalIOError);
    }

    read(contextDict, codeDict);

    const word sha1Str(sha1_.str());

    codeSha1Name_ = codeName_ + '_' + sha1Str;

    codeDirName_ =
    (
        codeDirName.empty()
      ? word('_' + sha1Str)
      : codeDirName
    );

    varSubstitutions_.set("typeName", codeName_);
    varSubstitutions_.set("uniqueFunctionName", codeSha1Name_);
    varSubstitutions_.set("SHA1sum", sha1Str);
}


Foam::dynamicCode::dynamicCode
(
    const dictionary& contextDict,
    const word& codeName,
    const word& codeDirName,
    const wordList& codeKeys,
    const wordList& codeDictVars,
    const word& codeOptionsFileName,
    const wordList& compileFiles,
    const wordList& copyFiles
)
:
    dynamicCode
    (
        contextDict,
        contextDict,
        codeName,
        codeDirName,
        codeKeys,
        codeDictVars,
        codeOptionsFileName,
        compileFiles,
        copyFiles
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicCode::read
(
    const dictionary& contextDict,
    const dictionary& codeDict
)
{
    // Expand all dictionary entries. Note that this removes any leading or
    // trailing whitespace, which is necessary for compilation options, and
    // doesn't hurt for everything else
    List<const entry*> codePtrs(codeKeys_.size(), nullptr);
    codeKeySubstitutions_.clear();
    forAll(codeKeys_, i)
    {
        const word& key = codeKeys_[i];
        codePtrs[i] = codeDict.lookupEntryPtr(key, false, false);
        if (codePtrs[i])
        {
            codeStrings_[i] = verbatimString(codePtrs[i]->stream());
            stringOps::inplaceExpandCodeString
            (
                codeStrings_[i],
                contextDict, // Lookup variables from the context dictionary
                codeDictVars_[i]
            );
            codeKeySubstitutions_.insert(key, codeStrings_[i]);
        }
        else
        {
            codeKeySubstitutions_.insert(key, "");
        }
    }

    // Code options
    const entry* optionsPtr =
        codeDict.lookupEntryPtr("codeOptions", false, false);
    if (optionsPtr)
    {
        optionsString_ = verbatimString(optionsPtr->stream());
        stringOps::inplaceExpandCodeString
        (
            optionsString_,
            contextDict,
            word::null
        );
        options_ = stringOps::trim(optionsString_);
    }

    // Code libs
    const entry* libsPtr = codeDict.lookupEntryPtr("codeLibs", false, false);
    if (libsPtr)
    {
        libsString_ = verbatimString(libsPtr->stream());
        stringOps::inplaceExpandCodeString
        (
            libsString_,
            contextDict,
            word::null
        );
        libs_ = stringOps::trim(libsString_);
    }

    // Calculate SHA1 digest from all entries
    OSHA1stream os;
    forAllConstIter(HashTable<string>, codeKeySubstitutions_, iter)
    {
        os << iter();
    }
    os << options_ << libs_;
    sha1_ = os.digest();

    // Add line directives after calculating SHA1
    forAll(codeKeys_, i)
    {
        if (codePtrs[i])
        {
            const word& key = codeKeys_[i];
            addLineDirective
            (
                codeKeySubstitutions_[key],
                codePtrs[i]->startLineNumber(),
                codeDict.name()
            );
        }
    }
}


Foam::word Foam::dynamicCode::libraryBaseName(const fileName& libPath)
{
    word libName(libPath.name(true));
    libName.erase(0, 3);    // Remove leading 'lib' from name
    return libName;
}


Foam::fileName Foam::dynamicCode::resolveTemplate
(
    const fileName& templateName
)
{
    return findConfigFile
    (
        templateName,
        codeTemplateDirName,
        "system"
    );
}


bool Foam::dynamicCode::copyOrCreateFiles(const bool verbose) const
{
    if (verbose)
    {
        Info<< "Creating new library in " << libRelPath() << endl;
    }

    HashTable<string> filterVars(varSubstitutions_);

    // Collect all the filter variables
    forAllConstIter(HashTable<string>, codeKeySubstitutions_, iter)
    {
        filterVars.set(iter.key(), iter());
    }

    const label nFiles =
        compileFiles_.size() + copyFiles_.size();

    DynamicList<fileName> resolvedFiles(nFiles);
    DynamicList<fileName> badFiles(nFiles);

    // Resolve template, or add to bad-files
    dynamicCode::resolveTemplates
    (
        compileFiles_,
        resolvedFiles,
        badFiles
    );
    dynamicCode::resolveTemplates
    (
        copyFiles_,
        resolvedFiles,
        badFiles
    );

    if (!badFiles.empty())
    {
        FatalErrorInFunction
            << "Could not find the code template(s): "
            << badFiles << nl
            << exit(FatalError);
    }

    // Create dir
    const fileName outputDir(codePath());

    // Create dir
    mkDir(outputDir);

    // Copy/filter files
    forAll(resolvedFiles, fileI)
    {
        const fileName& srcFile = resolvedFiles[fileI];
        const fileName dstFile(outputDir/srcFile.name());

        if (verbose)
        {
            Info << "    Copying " << srcFile << " to " << dstFile << endl;
        }

        IFstream is(srcFile);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Failed opening " << srcFile
                << exit(FatalError);
        }

        OFstream os(dstFile);
        if (!os.good())
        {
            FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
        }

        // Copy lines while expanding variables
        dynamicCode::copyAndFilter(is, os, filterVars);
    }


    // Create Make/files + Make/options
    createMakeFiles();
    createMakeOptions();

    writeDigest();

    return true;
}


bool Foam::dynamicCode::wmakeLibso() const
{
    Foam::string configCmd("cmake");
    configCmd+=" -S"+this->codePath();
    configCmd+=" -B"+this->codePath()/"Make";
    configCmd+=" -DCMAKE_C_COMPILER=$WM_CC";
    configCmd+=" -DCMAKE_CXX_COMPILER=$WM_CXX";
    configCmd+=" -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE";
    configCmd+=" -DCMAKE_PREFIX_PATH=$MIKENO_BINARY_INSTALL_PREFIX";
    configCmd+=" --no-warn-unused-cli";
#ifndef FULLDEBUG
    configCmd+=" --log-level=ERROR";
#endif
    Info<< "Invoking " << configCmd << endl;
    if (Foam::system(configCmd)) {
        return false;
    }

    Foam::string buildCmd("cmake");
    buildCmd+=" --build "+this->codePath()/"Make";

    Info<<"Invoking "<<buildCmd<<endl;
    if(Foam::system(buildCmd)) {
        return false;
    }

    return true;
}


bool Foam::dynamicCode::upToDate() const
{
    const fileName file(digestFile());

    if (!exists(file, false, true) || SHA1Digest(IFstream(file)()) != sha1_)
    {
        return false;
    }

    return true;
}


void* Foam::dynamicCode::loadLibrary(const fileName& libPath) const
{
    // Cached access to dl libs.
    // Guarantees clean up upon destruction of Time.
    if (libs.open(libPath, false))
    {
        return libs.findLibrary(libPath);
    }
    else
    {
        // Uncached opening of libPath. Do not complain if cannot be loaded
        return dlOpen(libPath, false);
    }
}


void Foam::dynamicCode::createLibrary
(
    const dictionary& dict,
    const bool masterOnlyRead
) const
{
    const bool create =
        Pstream::master()
     || (regIOobject::fileModificationSkew <= 0);   // Not NFS

    if (create)
    {
        // Write files for new library
        if (!upToDate())
        {
            if (!copyOrCreateFiles(true))
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Failed writing files for" << nl
                    << libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        if (!wmakeLibso())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Failed wmake " << libRelPath() << nl
                << exit(FatalIOError);
        }
    }

    // All processes must wait for compile to finish
    // Only block if not master only reading of a global dictionary
    if
    (
       !masterOnlyRead
     && regIOobject::fileModificationSkew > 0
    )
    {
        const fileName libPath = this->libPath();

        // Determine and communicate the master file size. Scattering
        // blocks the other processes until the master has finished
        // compiling.
        off_t masterSize = Pstream::master() ? fileSize(libPath) : -1;
        Pstream::scatter(masterSize);

        // Determine the local file size. This may be incorrect if NFS is
        // taking its time, in which case we wait and try again.
        off_t mySize = Pstream::master() ? masterSize : fileSize(libPath);

        if (debug)
        {
            Pout<< endl<< "on processor " << Pstream::myProcNo()
                << " have masterSize:" << masterSize
                << " and localSize:" << mySize
                << endl;
        }

        if (mySize < masterSize)
        {
            if (debug)
            {
                Pout<< "Local file " << libPath
                    << " not of same size (" << mySize
                    << ") as master ("
                    << masterSize << "). Waiting for "
                    << regIOobject::fileModificationSkew
                    << " seconds." << endl;
            }
            sleep(regIOobject::fileModificationSkew);

            // Recheck local size
            mySize = Foam::fileSize(libPath);

            if (mySize < masterSize)
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Cannot read (NFS mounted) library " << nl
                    << libPath << nl
                    << "on processor " << Pstream::myProcNo()
                    << " detected size " << mySize
                    << " whereas master size is " << masterSize
                    << " bytes." << nl
                    << "If your case is not NFS mounted"
                    << " (so distributed) set fileModificationSkew"
                    << " to 0"
                    << exit(FatalIOError);
            }
        }

        if (debug)
        {
            Pout<< endl<< "on processor " << Pstream::myProcNo()
                << " after waiting: have masterSize:" << masterSize
                << " and localSize:" << mySize
                << endl;
        }
    }
}


void Foam::dynamicCode::read(const dictionary& contextDict)
{
    read(contextDict, contextDict);
}


void Foam::dynamicCode::write(Ostream& os) const
{
    writeEntry(os, "name", codeName_);

    forAll(codeStrings_, i)
    {
        if (codeStrings_[i] != verbatimString::null)
        {
            writeEntry(os, codeKeys_[i], codeStrings_[i]);
        }
    }

    if (optionsString_ != verbatimString::null)
    {
        writeEntry(os, "codeOptions", optionsString_);
    }

    if (libsString_ != verbatimString::null)
    {
        writeEntry(os, "codeLibs", libsString_);
    }
}


// ************************************************************************* //
