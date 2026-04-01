#ifndef MIKENO_SRC_PARSE_WMAKE_H
#define MIKENO_SRC_PARSE_WMAKE_H

#include <string>
#include <string_view>
#include <map>
#include <vector>


namespace Foam::wmakeParse {

  [[nodiscard]] std::string_view trim(std::string_view str)noexcept;

  enum class undefined_reference_behavior {
    empty_string,
    throw_exception,
  };


  struct wmake_parse_option {
    undefined_reference_behavior when_undefined_reference=undefined_reference_behavior::empty_string;
    bool parse_brace_expression{false}; // Parse ${VAR}
  };

  std::vector<std::string> parse_wmake_file(std::string_view file_content,
                                            std::map<std::string,std::string>& parent_dict,
                                            const wmake_parse_option&option={});

  [[maybe_unused]] [[nodiscard]] std::map<std::string,std::string> get_environment_variables() noexcept;

  [[nodiscard]] std::string parse_link_libs(std::string_view compile_command) noexcept;

  [[nodiscard]] std::string parse_include_dirs(std::string_view raw_compile_command) noexcept;
}

#endif