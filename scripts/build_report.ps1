param([switch]$InstallMissingPackages)
$ErrorActionPreference = 'Stop'
$repositoryRoot = Split-Path -Parent $PSScriptRoot
$reportDirectory = Join-Path $repositoryRoot 'latex'
$outputDirectory = Join-Path $repositoryRoot 'output\pdf'
$pdfLaTeX = (Get-Command pdflatex -ErrorAction Stop).Source
$bibTeX = (Get-Command bibtex -ErrorAction Stop).Source
$installerOption = if ($InstallMissingPackages) { '-enable-installer' } else { '-disable-installer' }
# The installer switches are MiKTeX-specific; TeX Live needs neither.
$isMiKTeX = (& $pdfLaTeX --version | Select-Object -First 1) -match 'MiKTeX'
$latexArguments = @('-interaction=nonstopmode','-halt-on-error','-output-directory=build','main.tex')
if ($isMiKTeX) { $latexArguments = @($installerOption) + $latexArguments }
New-Item -ItemType Directory -Force -Path (Join-Path $reportDirectory 'build'),$outputDirectory | Out-Null
Push-Location $reportDirectory
try {
    & $pdfLaTeX @latexArguments
    if ($LASTEXITCODE -ne 0) { throw 'First LaTeX pass failed; see latex/build/main.log.' }
    & $bibTeX 'build/main'
    if ($LASTEXITCODE -ne 0) { throw 'BibTeX failed; see latex/build/main.blg.' }
    foreach ($pass in 1..2) {
        & $pdfLaTeX @latexArguments
        if ($LASTEXITCODE -ne 0) { throw "LaTeX pass $pass failed; see latex/build/main.log." }
    }
    $buildLog = Get-Content -LiteralPath 'build/main.log' -Raw
    if ($buildLog -match 'There were undefined (references|citations)|multiply-defined labels|Label\(s\) may have changed') {
        throw 'The PDF has unresolved cross-references; inspect latex/build/main.log.'
    }
    Copy-Item -LiteralPath 'build/main.pdf' -Destination (Join-Path $outputDirectory 'scipio1_design_report.pdf') -Force
    Write-Output "Report: $(Join-Path $outputDirectory 'scipio1_design_report.pdf')"
} finally {
    Pop-Location
}
