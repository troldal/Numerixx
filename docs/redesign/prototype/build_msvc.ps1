# usage: build_msvc.ps1 <source.cpp> <tag> [extra defines...]  -- MSVC cl 19.51 and clang-cl 22 (LLVM22), timed
param([string]$src, [string]$tag, [string[]]$extra = @())
. $PSScriptRoot\msvc_env.ps1
Set-Location $PSScriptRoot; New-Item -ItemType Directory -Force out | Out-Null
$eigen = if ($env:EIGEN_DIR) { $env:EIGEN_DIR } else { "..\eigen-5.0.1" }
$inc = @("/I.", "/IC:\Dev\XLThermo\FXT\include", "/I$eigen")
function Run-One($name, $exe, $args2) {
  $log = "out\${tag}_$name.log"
  $sw = [Diagnostics.Stopwatch]::StartNew()
  & $exe @args2 *> $log
  $rc = $LASTEXITCODE
  $sw.Stop()
  $warn = (Select-String -Path $log -Pattern "warning" -SimpleMatch | Measure-Object).Count
  $line = "{0,-22} compile rc={1} time={2,5:N2}s warnings={3}" -f $name, $rc, $sw.Elapsed.TotalSeconds, $warn
  if ($rc -eq 0) {
    $out = & ".\out\${tag}_$name.exe" 2>&1
    $rrc = $LASTEXITCODE
    $out | Set-Content "out\${tag}_$name.run"
    $line += "  run rc=$rrc  " + ($out | Select-Object -Last 1)
    $line
  } else {
    $line
    Select-String -Path $log -Pattern "error" | Select-Object -First 6 | ForEach-Object { "    " + $_.Line }
  }
}
$d = $extra | ForEach-Object { "/D$_" }
Run-One "msvc" "cl" (@("/nologo", "/std:c++latest", "/EHsc", "/permissive-", "/W4", "/O2", "/utf-8", "/Zc:__cplusplus") + $d + $inc + @($src, "/Fe:out\${tag}_msvc.exe", "/Fo:out\${tag}_msvc.obj"))
Run-One "clang-cl" $CLANGCL (@("--driver-mode=cl", "--no-default-config", "/std:c++latest", "/EHsc", "/W4", "/O2") + $d + $inc + @($src, "/Fe:out\${tag}_clang-cl.exe", "/Fo:out\${tag}_clang-cl.obj"))
