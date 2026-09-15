param(
  [string]$Remote = "git@github.com:fuyucheng514-tech/cellbit.git",
  [string]$TargetBranch = "microsags-v0.1.0"
)

$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

Write-Host "This will replace only the selected remote branch with this Microsags repository." -ForegroundColor Yellow
Write-Host "Remote: $Remote" -ForegroundColor Yellow
Write-Host "Target branch: $TargetBranch" -ForegroundColor Yellow
$answer = Read-Host "Type REPLACE to continue"
if ($answer -cne 'REPLACE') {
  Write-Host 'Cancelled.'
  exit 1
}

git remote remove origin 2>$null
git remote add origin $Remote
git push --force origin "main:$TargetBranch"
Write-Host "Published Microsags to branch $TargetBranch. The old remote main branch was not changed." -ForegroundColor Green
