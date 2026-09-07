param([string]$Python = 'python', [string]$OutputDirectory = '')
$ErrorActionPreference = 'Stop'
$releaseRoot = $PSScriptRoot
if (-not $OutputDirectory) {
    $OutputDirectory = Join-Path $releaseRoot ('outputs/replay_' + [guid]::NewGuid().ToString('N'))
}
if (Test-Path -LiteralPath $OutputDirectory) { throw 'Output directory exists; choose a new path.' }
& $Python -B (Join-Path $releaseRoot 'scripts/run_analysis.py') --out (Join-Path $OutputDirectory 'full827')
if ($LASTEXITCODE -ne 0) { throw 'Full baseline or feedback differs; author review required.' }
& $Python -B (Join-Path $releaseRoot 'scripts/run_analysis.py') --scope top50 --out (Join-Path $OutputDirectory 'top50')
if ($LASTEXITCODE -ne 0) { throw 'Top50 reproduction failed.' }
& $Python -B (Join-Path $releaseRoot 'scripts/figure6.py') --source (Join-Path $OutputDirectory 'full827') --out (Join-Path $OutputDirectory 'Figure6')
if ($LASTEXITCODE -ne 0) { throw 'Figure export failed.' }
& $Python -B -m pytest -p no:cacheprovider (Join-Path $releaseRoot 'tests') -q
if ($LASTEXITCODE -ne 0) { throw 'Release tests failed.' }
Write-Output ('Reproduction complete: ' + $OutputDirectory)
