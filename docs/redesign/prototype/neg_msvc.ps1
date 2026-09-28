. $PSScriptRoot\msvc_env.ps1
Set-Location $PSScriptRoot; New-Item -ItemType Directory -Force out | Out-Null
$tests = Get-ChildItem neg\neg_*.cpp
foreach ($t in $tests) {
  $b = $t.BaseName
  foreach ($c in @("msvc", "clang-cl")) {
    $log = "out\$b.$c.log"
    if ($c -eq "msvc") { cl /nologo /std:c++latest /EHsc /permissive- /W4 /Zs /I. /Ineg $t.FullName *> $log }
    else { & $CLANGCL --driver-mode=cl --no-default-config /std:c++latest /EHsc /W4 /Zs /I. /Ineg $t.FullName *> $log }
    $rc = $LASTEXITCODE
    $all = Get-Content $log
    $err = ($all | Select-String -Pattern "error" | Measure-Object).Count
    $first = ($all | Select-String -Pattern "error" | Select-Object -First 1).Line
    if ($first) { $first = ($first -replace '^.*?(error)', '$1'); if ($first.Length -gt 230) { $first = $first.Substring(0, 230) } }
    "{0,-26} {1,-8} rc={2} err-lines={3,-3} total-lines={4,-4} {5}" -f $b, $c, $rc, $err, $all.Count, $first
  }
}
"== delete-with-reason on cl:"
cl /nologo /std:c++latest /Zs neg\probe_delete_reason.cpp 2>&1 | Select-String "error|warning" | ForEach-Object { $_.Line -replace '^.*?(error|warning)', '$1' }
