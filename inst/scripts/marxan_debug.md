convert_eoln: system2(paste0("perl -pe 's/\r\n|\n|\r/\n/g' ",sInFile," > convert.dat"),wait=T)
[marxanui/ingest_marxan_data.R at master · mattwatts/marxanui](https://github.com/mattwatts/marxanui/blob/master/R/ingest_marxan_data.R#L144)

```
> launch_app("import")
hello
/Library/Frameworks/R.framework/Versions/3.5/Resources/library/marxanui/import

Listening on http://127.0.0.1:5789
marxanui_start start localuser
marxanui_start end
sUserIP
fingerprint:
ip:
userhostname: Invalid IP address
protocol: http:
hostname: 127.0.0.1
pathname: /
port: 5789
search:
queries:
puid_choices a puid_choices b
output$usermessages input$updateusermessages 0
sUserIP 128.111.61.49
fingerprint: 1991249ee558b321df87e7e97b1152c5
ip: 128.111.61.49
userhostname: 128-111-61-49.vpn.ucsb.edu
protocol: http:
hostname: 127.0.0.1
pathname: /
port: 5789
search:
queries:
output$usermessages input$updateusermessages 1
temp path /Users/bbest/marxanui/temp//upload_330151
unzip done
output$usermessages input$updateusermessages 2
Warning in system2(paste0("perl -pe 's/\r\n|\n|\r/\n/g' ", sInFile, " > convert.dat"),  :sh: line 3: perl -pe 's/
|
/
/g' /Users/bbest/marxanui/temp//upload_330151/scenario_01/input.dat > convert.dat: No such file or directory

  error in running command
sh: cp convert.dat /Users/bbest/marxanui/temp//upload_330151/scenario_01/input.dat: No such file or directory
Warning in system2(paste0("cp convert.dat ", sInFile), wait = T) :
  error in running command
reading input.dat /Users/bbest/marxanui/temp//upload_330151/scenario_01/input.dat
input.dat read /Users/bbest/marxanui/temp//upload_330151/scenario_01/input.dat
before fBoundDat
after fBoundDat
paths searched
sh: line 3: perl -pe 's/
|
/
/g' /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu.csv > convert.dat: No such file or directory
Warning in system2(paste0("perl -pe 's/\r\n|\n|\r/\n/g' ", sInFile, " > convert.dat"),  :
  error in running command
sh: cp convert.dat /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu.csv: No such file or directory
Warning in system2(paste0("cp convert.dat ", sInFile), wait = T) :
  error in running command
sh: line 3: perl -pe 's/
|
/
/g' /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/spp.csv > convert.dat: No such file or directory
Warning in system2(paste0("perl -pe 's/\r\n|\n|\r/\n/g' ", sInFile, " > convert.dat"),  :
  error in running command
sh: cp convert.dat /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/spp.csv: No such file or directory
Warning in system2(paste0("cp convert.dat ", sInFile), wait = T) :
  error in running command
sh: line 3: perl -pe 's/
|
/
/g' /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu_spp.csv > convert.dat: No such file or directory
Warning in system2(paste0("perl -pe 's/\r\n|\n|\r/\n/g' ", sInFile, " > convert.dat"),  :
  error in running command
sh: cp convert.dat /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu_spp.csv: No such file or directory
Warning in system2(paste0("cp convert.dat ", sInFile), wait = T) :
  error in running command
sh: line 3: perl -pe 's/
|
/
/g' /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/bound.csv > convert.dat: No such file or directory
Warning in system2(paste0("perl -pe 's/\r\n|\n|\r/\n/g' ", sInFile, " > convert.dat"),  :
  error in running command
sh: cp convert.dat /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/bound.csv: No such file or directory
Warning in system2(paste0("cp convert.dat ", sInFile), wait = T) :
  error in running command
reading input files
smart_read reading file /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu.csv
smart_read file read /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu.csv
smart_read reading file /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/spp.csv
smart_read file read /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/spp.csv
smart_read reading file /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu_spp.csv
smart_read file read /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/pu_spp.csv
smart_read reading file /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/bound.csv
smart_read file read /Users/bbest/marxanui/temp//upload_330151/scenario_01/input/bound.csv
marxan files ok
input files read
marxan files processed
pulayer dbf read
iupdateusermessages 3
OGR data source with driver: ESRI Shapefile
Source: "/Users/bbest/marxanui/temp/upload_330151/scenario_01/pulayer", layer: "pu"
with 113401 features
It has 1 fields
pulayer read /Users/bbest/marxanui/temp//upload_330151/scenario_01/pulayer/pu.shp
planning units simplified
planning units dissolved
```
3:00 am

...

```
pulayer.Rdata created /Users/bbest/marxanui/temp//upload_330151/marxan/pulayer/pulayer.Rdata
shapefiles processed
ParseMarxanZip end
length(ErrorMsg) 0
iupdateusermessages 4
output$usermessages input$updateusermessages 3
output$usermessages input$updateusermessages 4
acceptclicked
acceptclicked
click acceptupload 1
sDatabasePath /Users/bbest/marxanui//marxan/scenario_01
output$usermessages input$updateusermessages 5
```
