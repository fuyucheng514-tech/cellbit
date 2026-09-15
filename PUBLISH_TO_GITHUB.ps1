param(
  [string]$Remote = "git@github.com:fuyucheng514-tech/Microsags.git"
)

$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

Write-Host "This will replace the remote main branch contents with this Microsags repository." -ForegroundColor Yellow
Write-Host "Remote: $Remote" -ForegroundColor Yellow
$answer = Read-Host "Type REPLACE to continue"
if ($answer -cne 'REPLACE') {
  Write-Host 'Cancelled.'
  exit 1
}

git remote remove origin 2>$null
git remote add origin $Remote
git push --force origin main
Write-Host 'Published Microsags to GitHub.' -ForegroundColor Green
