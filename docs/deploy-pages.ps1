#!/usr/bin/env pwsh
<#
.SYNOPSIS
Build the Flexia manual and publish the compiled files to the "pages" branch.

.DESCRIPTION
1. Renders the Quarto book in this directory (HTML site, PDF, DOCX).
2. Clones the repository into a temporary directory, checks out (or creates)
   the "pages" branch, replaces its contents with the built files and pushes.

The working copy is never modified: the branch is prepared in a throw-away
clone under the system temporary directory.

.PARAMETER SkipRender
Publish the existing "docs/_book" directory without rendering it again.

.PARAMETER Remote
Name of the git remote that receives the "pages" branch (default: origin).

.PARAMETER Branch
Name of the publication branch (default: pages).

.PARAMETER Message
Commit message for the publication commit.

.PARAMETER NoPush
Prepare and commit the branch in the temporary clone but do not push and do
not delete the clone; the path to it is printed at the end.

.EXAMPLE
./docs/deploy-pages.ps1

.EXAMPLE
./docs/deploy-pages.ps1 -SkipRender -Message "docs: publish 1.0 manual"
#>
[CmdletBinding()]
param(
    [switch]$SkipRender,
    [string]$Remote = "origin",
    [string]$Branch = "pages",
    [string]$Message,
    [switch]$NoPush
)

$ErrorActionPreference = "Stop"
# Native git commands report failures through $LASTEXITCODE; do not turn them
# into terminating errors before we can inspect the exit code.
$PSNativeCommandUseErrorActionPreference = $false

$DocsDir = $PSScriptRoot
$BookDir = Join-Path $DocsDir "_book"

function Invoke-Git {
    param(
        [Parameter(Mandatory = $true)][string]$Path,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )

    & git -C $Path @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "git $($Arguments -join ' ') failed in '$Path' (exit code $LASTEXITCODE)"
    }
}

function Get-GitText {
    param(
        [Parameter(Mandatory = $true)][string]$Path,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )

    $value = & git -C $Path @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "git $($Arguments -join ' ') failed in '$Path' (exit code $LASTEXITCODE)"
    }
    return ($value | Out-String).Trim()
}

# --- проверки ---------------------------------------------------------------

$RepoRoot = Get-GitText -Path $DocsDir -Arguments @("rev-parse", "--show-toplevel")
if (-not $RepoRoot) {
    throw "'$DocsDir' is not inside a git repository"
}

$OriginUrl = Get-GitText -Path $RepoRoot -Arguments @("remote", "get-url", $Remote)

if (-not $SkipRender) {
    $quarto = Get-Command quarto -ErrorAction SilentlyContinue
    if (-not $quarto) {
        throw "quarto was not found in PATH; install Quarto or pass -SkipRender"
    }

    Write-Host "==> Rendering the manual (quarto render docs)"
    & quarto render $DocsDir
    if ($LASTEXITCODE -ne 0) {
        throw "quarto render failed (exit code $LASTEXITCODE)"
    }
}

if (-not (Test-Path -LiteralPath $BookDir)) {
    throw "Build directory '$BookDir' does not exist; run 'quarto render docs' first"
}

$required = @("index.html", "Flexia-manual.pdf", "Flexia-manual.docx")
foreach ($name in $required) {
    $path = Join-Path $BookDir $name
    if (-not (Test-Path -LiteralPath $path) -or (Get-Item -LiteralPath $path).Length -eq 0) {
        throw "Missing build artifact '$name' in '$BookDir'; run 'quarto render docs' first"
    }
}

# --- публикация -------------------------------------------------------------

$TempDir = Join-Path ([System.IO.Path]::GetTempPath()) ("flexia-pages-" + [guid]::NewGuid().ToString("N"))
$KeepTemp = $false

try {
    Write-Host "==> Cloning the repository to $TempDir"
    Invoke-Git -Path $DocsDir -Arguments @("clone", "--quiet", $RepoRoot, $TempDir)

    Invoke-Git -Path $TempDir -Arguments @("remote", "set-url", "origin", $OriginUrl)

    & git -C $TempDir fetch --quiet origin $Branch 2>$null
    if ($LASTEXITCODE -eq 0) {
        Invoke-Git -Path $TempDir -Arguments @("checkout", "--quiet", "-B", $Branch, "FETCH_HEAD")
    }
    else {
        Write-Host "==> Branch '$Branch' does not exist yet; creating an orphan branch"
        Invoke-Git -Path $TempDir -Arguments @("checkout", "--quiet", "--orphan", $Branch)
    }

    Write-Host "==> Replacing branch content with the compiled manual"
    Get-ChildItem -LiteralPath $TempDir -Force |
        Where-Object { $_.Name -ne ".git" } |
        ForEach-Object { Remove-Item -LiteralPath $_.FullName -Recurse -Force }

    Get-ChildItem -LiteralPath $BookDir -Force |
        ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination $TempDir -Recurse -Force }

    # GitHub Pages: не прогонять сайт через Jekyll.
    New-Item -ItemType File -Path (Join-Path $TempDir ".nojekyll") -Force | Out-Null

    # --force is required: the repository .gitignore ignores *.html and
    # docs/_book, but the publication branch must contain the built HTML.
    Invoke-Git -Path $TempDir -Arguments @("add", "--all", "--force")

    $changes = & git -C $TempDir status --porcelain
    if ($changes) {
        $userName = Get-GitText -Path $RepoRoot -Arguments @("config", "user.name")
        $userEmail = Get-GitText -Path $RepoRoot -Arguments @("config", "user.email")

        if (-not $userName -or -not $userEmail) {
            throw "git user.name and user.email are not configured; set them before publishing"
        }

        if (-not $Message) {
            $Message = "docs: publish compiled manual $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
        }

        Invoke-Git -Path $TempDir -Arguments @(
            "-c", "user.name=$userName",
            "-c", "user.email=$userEmail",
            "commit", "--quiet", "-m", $Message
        )
    }
    else {
        Write-Host "==> No changes in the compiled manual; nothing to commit"
    }

    if ($NoPush) {
        $KeepTemp = $true
        Write-Host "==> -NoPush: branch '$Branch' was prepared in $TempDir"
    }
    else {
        Write-Host "==> Pushing '$Branch' to '$Remote' ($OriginUrl)"
        Invoke-Git -Path $TempDir -Arguments @("push", "--quiet", "origin", "$Branch")
    }
}
finally {
    if (-not $KeepTemp -and (Test-Path -LiteralPath $TempDir)) {
        Remove-Item -LiteralPath $TempDir -Recurse -Force
    }
}

# --- итог -------------------------------------------------------------------

$pagesUrl = $null
if ($OriginUrl -match 'github\.com[:/](?<owner>[^/]+)/(?<repo>[^/]+?)(\.git)?$') {
    $pagesUrl = "https://$($Matches['owner']).github.io/$($Matches['repo'])/"
}

Write-Host ""
Write-Host "Done: '$Branch' was published to '$Remote'." -ForegroundColor Green
if ($pagesUrl) {
    Write-Host "GitHub Pages URL (if Pages is enabled for the '$Branch' branch): $pagesUrl"
}
