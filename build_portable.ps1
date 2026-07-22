[CmdletBinding()]
param(
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$buildRoot = Join-Path $root "build\portable"
$pyInstallerDist = Join-Path $root "dist\pyinstaller-portable"
$releaseRoot = Join-Path $root "dist\Ophiuchus_Portable"
$zipPath = Join-Path $root "dist\Ophiuchus_Portable.zip"
$healthPath = Join-Path $buildRoot "portable_health.json"

function Assert-ProjectChild([string]$PathToCheck) {
    $projectPrefix = [IO.Path]::GetFullPath($root).TrimEnd('\') + '\'
    $resolved = [IO.Path]::GetFullPath($PathToCheck)
    if (-not $resolved.StartsWith($projectPrefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to modify a path outside the project: $resolved"
    }
}

foreach ($path in @($buildRoot, $pyInstallerDist, $releaseRoot, $zipPath)) {
    Assert-ProjectChild $path
    if (Test-Path -LiteralPath $path) {
        Remove-Item -LiteralPath $path -Recurse -Force
    }
}
New-Item -ItemType Directory -Path $buildRoot, $pyInstallerDist, $releaseRoot | Out-Null

& $Python -m PyInstaller --version | Out-Host
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller is not installed in the selected build Python: $Python"
}

& $Python -m PyInstaller `
    --noconfirm `
    --clean `
    --workpath $buildRoot `
    --distpath $pyInstallerDist `
    (Join-Path $root "packaging\OphiuchusPortable.spec")
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller failed with exit code $LASTEXITCODE."
}

$frozenSource = Join-Path $pyInstallerDist "OphiuchusApp"
$runtimeTarget = Join-Path $releaseRoot "runtime"
if (-not (Test-Path -LiteralPath (Join-Path $frozenSource "OphiuchusApp.exe") -PathType Leaf)) {
    throw "Frozen OphiuchusApp.exe was not produced."
}
Copy-Item -LiteralPath $frozenSource -Destination $runtimeTarget -Recurse

& (Join-Path $root "build_launcher_exe.ps1")
$releaseFiles = @(
    "Ophiuchus.exe",
    "README.md",
    "VERSION",
    "install_desktop_shortcut.bat",
    "install_desktop_shortcut.ps1"
)
foreach ($relative in $releaseFiles) {
    Copy-Item -LiteralPath (Join-Path $root $relative) -Destination (Join-Path $releaseRoot $relative)
}
New-Item -ItemType Directory -Path (Join-Path $releaseRoot "docs") | Out-Null
Copy-Item -LiteralPath (Join-Path $root "docs\Ophiuchus_操作手册.md") -Destination (Join-Path $releaseRoot "docs\Ophiuchus_操作手册.md")

& (Join-Path $runtimeTarget "OphiuchusApp.exe") --health-check $healthPath
if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $healthPath -PathType Leaf)) {
    throw "Frozen runtime health check failed."
}
$health = Get-Content -LiteralPath $healthPath -Raw | ConvertFrom-Json
if (-not $health.core_ready -or -not $health.manual_available -or -not $health.user_data_writable) {
    throw "Frozen runtime health report did not pass the release gate."
}

$manifest = Join-Path $releaseRoot "SHA256SUMS.txt"
$hashLines = Get-ChildItem -LiteralPath $releaseRoot -Recurse -File |
    Where-Object { $_.FullName -ne $manifest } |
    Sort-Object FullName |
    ForEach-Object {
        $relative = $_.FullName.Substring($releaseRoot.Length).TrimStart('\').Replace('\', '/')
        $hash = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
        "$hash  $relative"
    }
[IO.File]::WriteAllLines($manifest, $hashLines, [Text.UTF8Encoding]::new($false))

Compress-Archive -LiteralPath $releaseRoot -DestinationPath $zipPath -CompressionLevel Optimal
$zipHash = (Get-FileHash -LiteralPath $zipPath -Algorithm SHA256).Hash.ToLowerInvariant()
[IO.File]::WriteAllText("$zipPath.sha256", "$zipHash  Ophiuchus_Portable.zip`r`n", [Text.UTF8Encoding]::new($false))

Write-Host "Portable release: $releaseRoot"
Write-Host "Portable archive: $zipPath"
Write-Host "Health report: $healthPath"
