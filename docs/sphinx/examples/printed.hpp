// What a documentation program prints, captured so that its main can compare it with the output its page quotes.
// capture(argv[0]) sends stdout to a file named after the program in the current directory (under node, argv[0] is a
// host path that Emscripten's file system does not have); printed() reads it back, and writes it to stderr too, where
// the test's log shows it. A page quotes the expected output from a raw string in its program, which starts on
// the line after R"( so that the page can include its lines: quoted() drops that first newline.

#pragma once

#include <cstdio>
#include <string>
#include <string_view>

namespace docs {

    inline bool capture(char const *program)
    {
        std::string const name = program;
        std::string const path = name.substr(name.find_last_of("/\\") + 1) + ".out";
#if defined(_MSC_VER)
        // The Microsoft runtime deprecates freopen in favour of its freopen_s.
        std::FILE *file = nullptr;
        return freopen_s(&file, path.c_str(), "w+", stdout) == 0;
#else
        return std::freopen(path.c_str(), "w+", stdout) != nullptr;
#endif
    }

    inline std::string printed()
    {
        std::fflush(stdout);
        std::rewind(stdout);
        std::string text;
        for (int c = std::getc(stdout); c != EOF; c = std::getc(stdout)) text += static_cast<char>(c);
        std::fputs(text.c_str(), stderr);
        return text;
    }

    constexpr std::string_view quoted(std::string_view raw) { return raw.substr(1); }

} // namespace docs
