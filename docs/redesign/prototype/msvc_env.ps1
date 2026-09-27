# Import the MSVC 14.51 x64 environment into the current PowerShell session.
$vcvars = "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
$lines = cmd /c "`"$vcvars`" -vcvars_ver=14.51 > nul && set"
foreach ($l in $lines) { if ($l -match '^([^=]+)=(.*)$') { Set-Item -Path "env:$($matches[1])" -Value $matches[2] } }
$global:CLANGCL = "C:\Toolchains\LLVM22\bin\clang-22.exe"
