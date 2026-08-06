$file = $args[0]
#$file = "/nfs/APL_Genomics/raw/260804_N_I_080"
$drive = $file -split "\\"
$drive = $drive[0]
$drives = net use
$drive = ($drives -match ".*" + $drive + ".*")
$parts = $drive -split "\s{1,11}"
$windowsDrive = $parts[2]
$uncDrive = $parts[3]
$file = $file.replace($windowsDrive, $uncDrive).replace("\\healthy.bewell.ca\Apps\APL_Genomics", "\nfs\APL_Genomics").replace("\\healthy.bewell.ca\Apps\APL\Genomics_DEV", "\nfs\Genomics_DEV").replace("\","/")
$ip = "10.106.109.185" # sch01
$user_at_host="${env:UserName}@$ip"

$cmd="/nfs/Genomics_DEV/projects/alindsay/Projects/ngs-pipeline-launcher/ngs_pipeline_launcher/pipelineLauncherFromWin.sh $file"
ssh -t "$USER_AT_HOST" "bash --rcfile ~/.bashrc -c '$cmd'"
pause