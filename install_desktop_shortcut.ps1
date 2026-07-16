[CmdletBinding()]
param(
    [string]$ShortcutName = "Start Ophiuchus"
)

$ErrorActionPreference = "Stop"

$exeLauncher = Join-Path $PSScriptRoot "Ophiuchus.exe"
$batchLauncher = Join-Path $PSScriptRoot "start_ophiuchus.bat"
if (Test-Path -LiteralPath $exeLauncher -PathType Leaf) {
    $launcher = $exeLauncher
} elseif (Test-Path -LiteralPath $batchLauncher -PathType Leaf) {
    $launcher = $batchLauncher
} else {
    throw "No Ophiuchus launcher was found beside this installer."
}

$shell = New-Object -ComObject WScript.Shell
$desktop = $shell.SpecialFolders.Item("Desktop")
if ([string]::IsNullOrWhiteSpace($desktop)) {
    throw "Windows desktop folder could not be resolved."
}

$shortcutPath = Join-Path $desktop ($ShortcutName + ".lnk")
$shortcut = $shell.CreateShortcut($shortcutPath)
$shortcut.TargetPath = $launcher
$shortcut.WorkingDirectory = $PSScriptRoot
$shortcut.Description = "Launch Ophiuchus from its current folder"
$shortcut.Save()

Write-Host "Created desktop shortcut: $shortcutPath"
Write-Host "Target: $launcher"
