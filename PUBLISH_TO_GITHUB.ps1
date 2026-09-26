param(
  [string]$Remote = "git@github.com:fuyucheng514-tech/cellbit.git",
  [string]$TargetBranch = "main"
)

$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

Write-Host "This will publish the current Microsags commit to the selected remote branch." -ForegroundColor Yellow
Write-Host "Remote: $Remote" -ForegroundColor Yellow
Write-Host "Target branch: $TargetBranch" -ForegroundColor Yellow
$answer = Read-Host "Type PUBLISH to continue"
if ($answer -cne 'PUBLISH') {
  Write-Host 'Cancelled.'
  exit 1
}

if (git remote get-url origin 2>$null) {
  git remote set-url origin $Remote
} else {
  git remote add origin $Remote
}
git push origin "HEAD:$TargetBranch"
Write-Host "Published Microsags to branch $TargetBranch without rewriting history." -ForegroundColor Green
