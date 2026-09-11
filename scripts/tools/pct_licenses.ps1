# pct_licenses.ps1 -- how many Parallel Computing Toolbox seats are in use.
#
#   .\scripts\tools\pct_licenses.ps1            # prints "in use / issued (free)"
#   .\scripts\tools\pct_licenses.ps1 -Quiet     # prints just the number in use
#   .\scripts\tools\pct_licenses.ps1 -Users     # also lists who holds them
#
# Queries the Mayo FlexLM server with lmutil (the network license manager
# download), which is the only reliable answer: MATLAB's own `license`
# function reports what THIS session has checked out, not what the server has
# left. Exit code is the number of free seats (0 = none), so a launcher can
# test it.
#
# Written 2026-09-11 when run_numerics_verification's parpool failed with
# "Maximum number of simultaneous users for this product reached" (15 of 15
# seats taken, none by us).

param(
    [switch]$Quiet,
    [switch]$Users,
    [string]$LmUtil = "$env:USERPROFILE\Downloads\mathworks_network_license_manager_win64\etc\win64\lmutil.exe",
    [string]$Server = "27000@rcf-lmgrd4.mayo.edu",
    [string]$Feature = "Distrib_Computing_Toolbox"
)

if (-not (Test-Path $LmUtil)) {
    Write-Error "lmutil not found at $LmUtil"
    exit 255
}

$out = & $LmUtil lmstat -f $Feature -c $Server 2>&1
$line = $out | Select-String -Pattern "Total of (\d+) licenses issued;\s+Total of (\d+) licenses in use"
if (-not $line) {
    Write-Error "Could not parse lmstat output:`n$($out -join "`n")"
    exit 254
}
$issued = [int]$line.Matches[0].Groups[1].Value
$inUse  = [int]$line.Matches[0].Groups[2].Value
$free   = $issued - $inUse

if ($Quiet) {
    Write-Output $inUse
} else {
    Write-Output ("{0} in use / {1} issued ({2} free)" -f $inUse, $issued, $free)
    if ($Users) {
        $out | Select-String -Pattern "^\s+\S+ \S+ \S+ \(v\d+\)" | ForEach-Object {
            $f = ($_.Line.Trim() -split "\s+")
            Write-Output ("  {0,-10} {1,-22} since {2}" -f $f[0], $f[1], ($_.Line -replace ".*start ", ""))
        }
    }
}
exit $free
