. $PSScriptRoot\msvc_env.ps1 2>$null
Set-Location $PSScriptRoot; New-Item -ItemType Directory -Force out | Out-Null
$eigen = if ($env:EIGEN_DIR) { $env:EIGEN_DIR } else { "..\eigen-5.0.1" }
$inc = @("/I.", "/IC:\Dev\XLThermo\FXT\include", "/I$eigen")
$tus = [ordered]@{ "base" = @("base.cpp"); "core_nopipes" = @("test_core.cpp", "/DNXX_NO_PIPES", "/DNXX_FIX_UNWRAP_FAILURE"); "core_pipes" = @("test_core.cpp", "/DNXX_FIX_UNWRAP_FAILURE"); "nd_inhouse" = @("test_nd.cpp"); "nd_eigen" = @("test_nd.cpp", "/DNXX_WITH_EIGEN"); "eigen_include_only" = @("eigen_only.cpp") }
function Best($exe, $a) { $m = 999.0; for ($i = 0; $i -lt 3; $i++) { $sw = [Diagnostics.Stopwatch]::StartNew(); & $exe @a *> $null; if ($LASTEXITCODE -ne 0) { return "FAIL" }; $sw.Stop(); $m = [Math]::Min($m, $sw.Elapsed.TotalSeconds) }; return "{0:N2}" -f $m }
"{0,-20} {1,8} {2,8} {3,8} {4,8}" -f "TU", "cl-O2", "cl-Od", "ccl-O2", "ccl-Od"
foreach ($k in $tus.Keys) {
  $s = $tus[$k]
  $a = Best "cl" (@("/nologo", "/std:c++latest", "/EHsc", "/permissive-", "/O2", "/c", "/Foout\tt.obj") + $inc + $s)
  $b = Best "cl" (@("/nologo", "/std:c++latest", "/EHsc", "/permissive-", "/Od", "/c", "/Foout\tt.obj") + $inc + $s)
  $c = Best $CLANGCL (@("--driver-mode=cl", "--no-default-config", "/std:c++latest", "/EHsc", "/O2", "/c", "/Foout\tt.obj") + $inc + $s)
  $d = Best $CLANGCL (@("--driver-mode=cl", "--no-default-config", "/std:c++latest", "/EHsc", "/Od", "/c", "/Foout\tt.obj") + $inc + $s)
  "{0,-20} {1,8} {2,8} {3,8} {4,8}" -f $k, $a, $b, $c, $d
}
