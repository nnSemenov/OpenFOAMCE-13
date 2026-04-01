//
// Created by joseph on 2025/11/15.
//

#include <ranges>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <optional>
#include <exception>
#include <cstdlib>
#include <format>

#include <unistd.h>

#include "parse_wmake.H"

namespace Foam::wmakeParse {

  std::string_view trim(std::string_view sv) noexcept {
    std::string_view empty=" \t\n\r\f\v";
    const size_t idx_beg=sv.find_first_not_of(empty);
    const size_t idx_end=sv.find_last_not_of(empty);
    if(idx_beg==std::string_view::npos) {
      assert(idx_end==std::string_view::npos);
      return {};
    }
    if(idx_end<idx_beg) {
      return {};
    }
    return {sv.data()+idx_beg,sv.data()+idx_end+1};
  }

  std::vector<std::string> filter_file_lines(std::string_view text) noexcept {
    std::vector<std::string> ret;
    std::string new_line;
    for(const auto & seg: text | std::views::split('\n')) {
      std::string_view raw_line{seg};
      std::string_view line= trim(raw_line);
      const bool is_empty_line = line.empty() or line.starts_with('#');

      if(line.ends_with('\\')) {
        if(not is_empty_line)
          new_line+=line.substr(0,line.length()-1);
        continue;
      }

      if(not is_empty_line)
        new_line+=line;

      auto trim_new_line=trim(new_line);
      if(not trim_new_line.empty()) {
        ret.emplace_back(trim_new_line);
      }
      new_line.clear();
    }


    auto trim_new_line=trim(new_line);
    if(not trim_new_line.empty()) {
      ret.emplace_back(trim_new_line);
    }
    return ret;
  }

  std::optional<size_t> first_quote_location(std::string_view line) {
    size_t quote_counter=0;
    for(size_t idx=0;idx<line.size();idx++) {
      const char ch=line[idx];
      if(ch=='\'' or ch=='\"') {
        quote_counter++;
      }
      if(ch=='=') {
        if(quote_counter>0) { // = after quote, not a variable statement
          return std::nullopt;
        }
        return idx;
      }
    }
    return std::nullopt;
  }

  std::string evaluate_string(std::string_view raw,
                              const std::map<std::string,std::string>& vars,
                              const wmake_parse_option&option) {
    std::string ret;
    for(const auto & seg_: raw | std::views::split('$')) {
      std::string_view seg{seg_};
      if(seg.empty()) {
        continue;
//        throw std::runtime_error{std::format("Found invalid expression: \"{}\"",raw)};
      }

      char bracket=0;
      char reverse_bracket=0;
      if(seg[0]=='(') {
        bracket=seg[0];
        reverse_bracket=')';
      }
      if(option.parse_brace_expression and seg[0]=='{') {
        bracket=seg[0];
        reverse_bracket='}';
      }

      if(bracket==0) {
        const ptrdiff_t head_offset=seg.data()-raw.data();
        if(head_offset<=0) { // Head of string ,somthing link EXE_INC = $(...)
          ret+=seg;
        }else { // Middle of string, like EXE_INC = $ENV_VAR
          ret+='$';
          ret+=seg;
        }
        continue;
      }

      const size_t reverse_bracket_loc=seg.find_first_of(reverse_bracket);
      if(reverse_bracket_loc==std::string_view::npos) {
        throw std::runtime_error{std::format("Found incomplete evaluate expression: \"${}\"",raw)};
      }

      const std::string var_name{seg.data()+1,seg.data()+reverse_bracket_loc};
      auto it = vars.find(var_name);
      const char* value = nullptr;
      if(it==vars.end()) { // non-existence variable
        switch (option.when_undefined_reference) {
          case undefined_reference_behavior::empty_string:
            value="";
            break;
          case undefined_reference_behavior::throw_exception:
            throw std::runtime_error{std::format("Undefined reference to variable \"{}\"",var_name)};
        }
      } else{
        value=it->second.c_str();
      }
      assert(value not_eq nullptr);

      ret+=value;
      std::string_view reset_part{seg.begin()+reverse_bracket_loc+1,seg.end()};
      ret += reset_part;
    }

    return ret;
  }

  std::vector<std::string> parse_wmake_file(std::string_view text, std::map<std::string,std::string>& parent_dict,
                                            const wmake_parse_option&option) {
    std::vector<std::string> files;
    auto lines = filter_file_lines(text);
    for(std::string & line: lines) {
      if(line.contains('$')) {
        line= evaluate_string(line,parent_dict,option);
      }

      const auto quote_loc_opt= first_quote_location(line);
      if(not quote_loc_opt) { // not a variable statement
        files.emplace_back(line);
        continue;
      }

      const size_t quote_loc=quote_loc_opt.value();
      std::string_view var_name_raw{line.data(),quote_loc};
      std::string_view var_name = trim(var_name_raw);

      std::string_view var_value_raw{line.begin()+ptrdiff_t(quote_loc+1),line.end()};
      std::string_view var_value= trim(var_value_raw);

      parent_dict.emplace(var_name,var_value);
    }
    return files;
  }

  [[maybe_unused]] std::map<std::string,std::string> get_environment_variables() noexcept {
    std::map<std::string,std::string> vars;

    for(size_t idx=0;;idx++) {
      const char* expression_c_str=environ[idx];
      if(expression_c_str== nullptr) {
        break;
      }
      std::string_view expression{expression_c_str};
      if(expression.empty()) {
        break;
      }
      const size_t eq_loc=expression.find_first_of('=');
      std::string var_name, value;
      if(eq_loc==std::string_view::npos) {
        value="";
        var_name=expression;
      }else {
        var_name=expression.substr(0,eq_loc);
        value=expression.substr(eq_loc+1,std::string_view ::npos);
      }

      vars.emplace(var_name,value);
    }

    return vars;
  }

  std::string parse_link_libs(std::string_view raw_compile_command) noexcept {
    std::string_view command=trim(raw_compile_command);
    std::string_view prefix="-l";
    if (not command.starts_with(prefix)) {
      return {};
    }
    return std::string{command.substr(prefix.size())};
  }

  std::string parse_include_dirs(std::string_view raw_compile_command) noexcept {
    std::string_view command=trim(raw_compile_command);
    for (std::string_view prefix :{"-I","-isystem","-i"}) {
      if (command.starts_with(prefix)) {
        return std::string{command.substr(prefix.size())};
      }
    }
    return {};
  }
}