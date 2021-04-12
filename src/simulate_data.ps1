[CmdletBinding()]
param (
    [Parameter(Mandatory)]
    [string]
    $omToolsJar,
    [Parameter(Mandatory)]
    [string]
    $ref,
    [switch]
    $noError
)
$ErrorActionPreference = "Stop"
$fileName = (Get-Item $ref).DirectoryName + "\" + (Get-Item $ref).BaseName + $(If ($noError) { "_without_error" } Else { "_with_error" })
$simulatedMap = $fileName + "_simulated.sdata"
# $simulatedMapCmap = $fileName + "_simulated.cmap"
$alignment = $fileName + "_alignment.omd"
$alignmentStats = $fileName + "_alignment_stats.txt"

$aptMapDataGeneratorArgs = @("--refmapin", $ref, "--optmapout", $simulatedMap, "--seed", "123")
if ($noError) {
    $aptMapDataGeneratorArgs += "--rsln"
    $aptMapDataGeneratorArgs += "0"
    $aptMapDataGeneratorArgs += "--meas"
    $aptMapDataGeneratorArgs += "0"
    $aptMapDataGeneratorArgs += "--scalesd"
    $aptMapDataGeneratorArgs += "0"
    $aptMapDataGeneratorArgs += "--subound"
    $aptMapDataGeneratorArgs += "1"
    $aptMapDataGeneratorArgs += "--slbound"
    $aptMapDataGeneratorArgs += "1"
    $aptMapDataGeneratorArgs += "--fpr"
    $aptMapDataGeneratorArgs += "0"
    $aptMapDataGeneratorArgs += "--fnr"
    $aptMapDataGeneratorArgs += "0"
}
# java -jar OMTools.jar FastaToOM --fastain ../ecoli.fasta --enzyme BspQI --refmapout ../ecoli.cmap
java -jar $omToolsJar OptMapDataGenerator $aptMapDataGeneratorArgs
java -jar $omToolsJar OMBlastMapper --refmapin $ref --optmapin $simulatedMap --optresout $alignment
((Get-Content -path $alignment -Raw) -replace ',', '.') | Set-Content -Path $alignment
java -jar $omToolsJar ResultStatistics --refmapin $ref --optmapin $simulatedMap --optresin $alignment --statout $alignmentStats
# java -jar $omToolsJar DataTools --optmapin $simulatedMap --optmapout $simulatedMapCmap
# ((Get-Content -path $simulatedMapCmap -Raw) -replace ',', '.') | Set-Content -Path $simulatedMapCmap