[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"

$source = Join-Path $PSScriptRoot "tools\OphiuchusLauncher.cs"
$output = Join-Path $PSScriptRoot "Ophiuchus.exe"
$compilerCandidates = @(
    (Join-Path $env:WINDIR "Microsoft.NET\Framework64\v4.0.30319\csc.exe"),
    (Join-Path $env:WINDIR "Microsoft.NET\Framework\v4.0.30319\csc.exe")
)
$compiler = $compilerCandidates | Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1

if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
    throw "Launcher source not found: $source"
}
if (-not $compiler) {
    throw "The Windows .NET Framework C# compiler was not found."
}
if (Test-Path -LiteralPath $output -PathType Leaf) {
    Remove-Item -LiteralPath $output -Force
}

$compilerArgs = @(
    "/nologo",
    "/target:winexe",
    "/optimize+",
    "/reference:System.Windows.Forms.dll",
    "/out:$output",
    $source
)
& $compiler @compilerArgs
if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $output -PathType Leaf)) {
    throw "Ophiuchus.exe build failed with exit code $LASTEXITCODE."
}

Write-Host "Built portable launcher: $output"
