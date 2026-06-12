/*
License
    This file is part of Mikeno. Mikeno (https://github.com/nnSemenov/Mikeno) is
    a fork of OpenFOAM (https://github.com/OpenFOAM/OpenFOAM-dev).

    Mikeno is free software: you can redistribute it and/or modify it  under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    Mikeno is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License
    along with Mikeno.  If not, see <http://www.gnu.org/licenses/>.

Application
    autoCompress

Description
    Compress all heavy fields (non-uniform, including mesh) to gzip, regardless
    of format(ascii or binary)
\*---------------------------------------------------------------------------*/
#include <zlib.h>
#include <filesystem>
#include <vector>
#include <set>
#include <format>
#include <cstdio>
#include <expected>

#include <omp.h>
#include <mutex>
#include <atomic>
#include <thread>

#include "Time.H"
#include "messageStream.H"
#include "IFstream.H"
#include "argList.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "fieldDictionary.H"

namespace stdfs = std::filesystem;

// Foam::word extract_file_type()
using namespace Foam;

//- Read dictionary from file and return
//  Sets stream to binary mode if specified in the optional header
struct file_info
{
    IOstream::streamFormat format;
    word file_class;
};

std::expected<file_info, std::string> get_dict_info(
    const fileName &dictFileName) noexcept;

struct foam_glob_option
{
    stdfs::directory_options fs_option{stdfs::directory_options::none};
    bool include_constant{false};
    bool include_zero{false};
};

std::vector<stdfs::path> glob_compressible_files(
    const Time &runTime, messageStream &info_dest,
    const std::set<std::string> &compressible_class,
    const foam_glob_option &option) noexcept;

struct compress_stat
{
    size_t original_bytes;
    size_t compressed_bytes;
};

std::expected<compress_stat, std::string> compress_file(
    const stdfs::path &file, const stdfs::path &dest, int level, bool keep_old,
    std::vector<uint8_t> &buf) noexcept;

int main(int argc, char **argv)
{
    using namespace Foam;
    argList::addOption("level", "0~9", "Compress level");
    argList::addOption("j", "0", "Number of threads for compilation");
    argList::addBoolOption("follow-symlink",
                           "Follow symbolic link when processing data files");
    argList::addBoolOption("list-class",
                           "Print all compressible classes and exit");
    argList::addBoolOption("withZero",
                           "include the '0/' dir in the times list");
    argList::addBoolOption("constant",
                           "include the 'constant/' dir in the times list");
    argList::addBoolOption("keep", "Keep original files");
    argList::noParallel();

#include "setRootCase.H"
#include "createTime.H"

    const bool follow_symlink = args.optionFound("follow-symlink");
    const auto fs_option = follow_symlink
        ? stdfs::directory_options::follow_directory_symlink
        : stdfs::directory_options::none;

    const int compress_level = args.optionLookupOrDefault<int>("level", 9);
    if (compress_level < 0 or compress_level > 9) {
        FatalErrorIn(args.executable())
            << std::format("invalid compress level for gzip: {}. Expected "
                           "integer in [0,9]",
                           compress_level)
                   .c_str()
            << exit(FatalError);
        return 1;
    }

    std::set<std::string> compressible_class{
        volScalarField::typeName,
        volVectorField::typeName,
        volSphericalTensorField::typeName,
        volSymmTensorField::typeName,
        volTensorField::typeName,

        surfaceScalarField::typeName,
        surfaceVectorField::typeName,
        surfaceSphericalTensorField::typeName,
        surfaceSymmTensorField::typeName,
        surfaceTensorField::typeName,

        pointScalarField::typeName,
        pointVectorField::typeName,
        pointSphericalTensorField::typeName,
        pointSymmTensorField::typeName,
        pointTensorField::typeName,

        volScalarField::Internal::typeName,
        volVectorField::Internal::typeName,
        volSphericalTensorField::Internal::typeName,
        volSymmTensorField::Internal::typeName,
        volTensorField::Internal::typeName,

        "scalarField",
        "vectorField",
        "labelList",
        "faceList",
        "regIOobject",
        "faceLabels",
        "flipMap",
        "cellLabels",
        "pointLabels",
        "faceCompactList",
        "cellZoneList",

    };

    {
        const bool print = args.optionFound("list-class");
        if (print) {
            for (const auto &class_ : compressible_class) {
                Info << class_.c_str() << '\n';
            }
            return 0;
        }
    }
    foam_glob_option glob_option;
    glob_option.fs_option = fs_option;
    glob_option.include_zero = args.optionFound("withZero");
    glob_option.include_constant = args.optionFound("constant");

    const auto compressible_files =
        glob_compressible_files(runTime, Info, compressible_class, glob_option);

    const bool compress_keep = args.optionFound("keep");

    const stdfs::path case_path = runTime.path().c_str();
    int num_threads = args.optionLookupOrDefault<int>("j", 1);
    if (num_threads <= 0) {
        num_threads = static_cast<int>(std::thread::hardware_concurrency());
    }
    Info << "Compressing with " << num_threads << " threads" << endl;
    omp_set_num_threads(num_threads);

    // Parallel compression
    std::mutex info_lock;
    std::atomic_int error_num{0};
#pragma omp parallel for schedule(runtime)
    // default(none)
    //  shared(compressible_files,compress_level,error_num,case_path,info_lock,Info,args,FatalError,compress_keep)
    for (const stdfs::path &src_file : compressible_files) {
        thread_local std::vector<uint8_t> buffer;
        if (error_num > 0) {
            continue;
        }

        stdfs::path compressed = src_file;
        compressed += ".gz";
        const auto relative_path = stdfs::relative(src_file, case_path);

        auto compress_opt = compress_file(src_file, compressed, compress_level,
                                          compress_keep, buffer);
        if (not compress_opt) {
            error_num++;
            FatalErrorIn(args.executable())
                << std::format("Compressing {} ... Failure: {}",
                               relative_path.string(), compress_opt.error())
                << endl;
            continue;
        }
        auto &compress_info = compress_opt.value();
        const scalar ratio = 100 * scalar(compress_info.compressed_bytes) /
            scalar(compress_info.original_bytes);

        {
            std::lock_guard<std::mutex> lkgd{info_lock};
            Info << std::format("Compressing {} ... Finish, {} bytes => {} "
                                "bytes, {:.2f}%",
                                relative_path.string(),
                                compress_info.original_bytes,
                                compress_info.compressed_bytes, ratio)
                        .c_str()
                 << endl;
        }
    }

    if (error_num > 0) {
        Info << std::format(
            "{} file(s) failed to be compressed. Exit with error",
            error_num.load());
        return 1;
    }

    return 0;
}

std::expected<file_info, std::string> get_dict_info(
    const fileName &dictFileName) noexcept
{
    IOstream::streamFormat dictFormat = IOstream::ASCII;

    file_info ret;
    // Read the first entry and if it is FoamFile set the file format

    IFstream dictFile(dictFileName);
    if (not dictFile().good()) {
        return std::unexpected{
            std::format("Can not open file {}", dictFileName.c_str())};
    }

    // Check if the first token in the file is "FoamFile"
    // to avoid problems if the first entry is a variable or function
    token firstToken;
    dictFile.read(firstToken);
    if ((not firstToken.isWord()) or
        (firstToken.wordToken() not_eq IOobject::foamFile))
        return std::unexpected{"not a foam file"};

    dictFile.putBack(firstToken);

    // Read the first entry from the dictionary
    autoPtr<entry> firstEntry(entry::New(dictFile()));
    if (not firstEntry->isDict()) {
        return std::unexpected{"FoamFile is not a dictionary"};
    }

    if (firstEntry->keyword() not_eq IOobject::foamFile) {
        return std::unexpected{"First entry is not the \"FoamFile\" header"};
    }

    const auto &header = firstEntry->dict();
    //        ret.file_class = firstEntry->dict().lookup<word>("class");
    if (not header.found("class")) {
        return std::unexpected{"FoamFile header doesn't specify class"};
    }
    ret.file_class = header.lookup<word>("class");

    if (not header.found("format")) {
        return std::unexpected{"FoamFile header doesn't specify format"};
    }
    const word format_word = header.lookup<word>("format");
    ret.format = IOstream::formatEnum(format_word);

    ret.format = dictFormat;
    return ret;
}

std::vector<stdfs::path> glob_compressible_files(
    const Time &runTime, messageStream &info_dest,
    const std::set<std::string> &compressible_class,
    const foam_glob_option &option) noexcept
{
    std::vector<stdfs::path> ret;
    const stdfs::path case_path = runTime.path().c_str();

    const auto times = runTime.times();
    for (const auto &instant : times) {
        if (instant.name() == "constant" and (not option.include_constant)) {
            continue;
        }
        if (instant.name() == "0" and (not option.include_zero)) {
            continue;
        }

        const word dir = runTime.path() / instant.name();
        //        Info<<std::format("Files under time =
        //        {}",instant.value()).c_str()<<endl;
        for (auto &entry : stdfs::recursive_directory_iterator{
                 dir.c_str(), option.fs_option}) {
            if (entry.is_directory()) {
                continue;
            }
            if (entry.is_symlink()) {
                info_dest << "\tSkip non-directory symlink " << entry.path()
                          << endl;
                continue;
            }
            if (entry.path().extension() == ".gz") {
                continue;
            }
            if (entry.path().extension() == ".orig") {
                continue;
            }
            const auto relative_path = relative(entry.path(), case_path);

            auto info_opt = get_dict_info(fileName{entry.path()});
            if (not info_opt) {
                info_dest << "Skip non-foam-file " << relative_path.c_str()
                          << endl;
                continue;
            }

            const file_info &info = info_opt.value();

            if (not compressible_class.contains(info.file_class)) {
                continue;
            }

            ret.emplace_back(entry.path());
        }
    }
    return ret;
}

std::expected<compress_stat, std::string> compress_file(
    const stdfs::path &file, const stdfs::path &dest, int level, bool keep_old,
    std::vector<uint8_t> &buf) noexcept
{
    FILE *ifile = fopen64(file.c_str(), "rb+");
    if (ifile == nullptr) {
        return std::unexpected{
            std::format("Failed to read file {}", file.string())};
    }

    level = std::clamp<int>(level, 0, 9);
    const std::string gz_mode = std::format("wb{}", level);

    gzFile ofile = gzopen64(dest.c_str(), gz_mode.data());
    if (ofile == nullptr) {
        return std::unexpected{
            std::format("Failed to write to {}", dest.string())};
    }

    const size_t N_buf = 8192;
    int gz_status = Z_OK;
    gz_status = gzbuffer(ofile, N_buf);
    if (gz_status not_eq Z_OK) {
        return std::unexpected{std::format(
            "Failed to set gzbuffer with error code {}", gz_status)};
    }

    buf.resize(N_buf);
    memset(buf.data(), 0, N_buf);
    size_t n_total_read = 0;
    //  size_t n_total_write=0;
    while (true) {
        const size_t n_read = fread(buf.data(), 1, N_buf, ifile);
        n_total_read += n_read;
        gzwrite(ofile, buf.data(), n_read);

        if (feof(ifile)) {
            break;
        }
    }

    fclose(ifile);
    gzclose(ofile);

    const size_t n_total_write = file_size(dest);

    if (not keep_old) {
        if (not stdfs::remove(file)) {
            return std::unexpected{
                std::format("Failed to remove old file {}", file.string())};
        }
    }

    return compress_stat{
        .original_bytes = n_total_read,
        .compressed_bytes = n_total_write,
    };
}