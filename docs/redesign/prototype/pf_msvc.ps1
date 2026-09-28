. $PSScriptRoot\msvc_env.ps1
Set-Location $PSScriptRoot
"== cl"; cl /nologo /std:c++latest /EHsc /permissive- probe_features.cpp /Fe:out\pf_msvc.exe /Fo:out\ | Out-Null; .\out\pf_msvc.exe
"== clang-cl 22"; & $CLANGCL --driver-mode=cl --no-default-config /std:c++latest /EHsc probe_features.cpp /Fe:out\pf_clangcl.exe /Fo:out\ ; .\out\pf_clangcl.exe
