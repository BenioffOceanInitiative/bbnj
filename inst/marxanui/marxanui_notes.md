
## marxanui

```
devtools::install_github("mattwatts/marxanui")
```

### attempt 1 error: `could not find function "marxanui_start"`

```
> Sys.getenv("HOME")
[1] "/Users/bbest"
> shiny::runApp('inst/marxan')
Loading required package: shiny
hello
/Users/bbest/github/marxanui/inst/marxan

Listening on http://127.0.0.1:3989
Loading required package: foreach
Loading required package: iterators
Loading required package: parallel
Loading required package: foreign
Loading required package: labdsv
Loading required package: mgcv
Loading required package: nlme
This is mgcv 1.8-25. For overview type 'help("mgcv-package")'.
Loading required package: MASS
Loading required package: cluster

Attaching package: ‘labdsv’

The following object is masked from ‘package:stats’:

    density

Loading required package: maptools
Loading required package: sp
Checking rgeos availability: TRUE
Loading required package: PBSmapping

-----------------------------------------------------------
PBS Mapping 2.70.5 -- Copyright (C) 2003-2018 Fisheries and Oceans Canada

PBS Mapping comes with ABSOLUTELY NO WARRANTY;
for details see the file COPYING.
This is free software, and you are welcome to redistribute
it under certain conditions, as outlined in the above file.

A complete user guide 'PBSmapping-UG.pdf' is located at
/Library/Frameworks/R.framework/Versions/3.5/Resources/library/PBSmapping/doc/PBSmapping-UG.pdf

Packaged on 2018-06-05
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
https://github.com/pbs-software

To see demos, type '.PBSfigs()'.
-----------------------------------------------------------


Loading required package: sqldf
Loading required package: gsubfn
Loading required package: proto
Loading required package: RSQLite
Loading required package: vegan
Loading required package: permute
Loading required package: lattice
This is vegan 2.5-3
Loading required package: xtable

Attaching package: ‘xtable’

The following object is masked from ‘package:maptools’:

    label

Warning: Error in marxanui_start: could not find function "marxanui_start"
  47: server [/Users/bbest/github/marxanui/inst/marxan/server.R#190]
Error in marxanui_start("marxan") :
  could not find function "marxanui_start"
sUserIP
sFingerprint
```

## working 2: after `library(marxanui)`

```
> runApp('inst/marxan')
hello
/Users/bbest/github/marxanui/inst/marxan

Listening on http://127.0.0.1:3989
Loading required package: marxanui

Attaching package: ‘marxanui’

The following object is masked from ‘package:base’:

    list.dirs

marxanui_start start localuser
create user home
trying URL 'http://marxan.net/downloads/unix_data.zip'
Content type 'application/zip' length 2856051 bytes (2.7 MB)
==================================================
downloaded 2.7 MB

marxanui_start end
sUserIP
sFingerprint
m
blm
input$zoomtoextent
input$zoomtoprev
fingerprint:
ip:
userhostname: Invalid IP address
protocol: http:
hostname: 127.0.0.1
pathname: /
port: 3989
search:
queries:

session$clientData
generate_aspect_width iAspectWidth 600
generate_aspect_height iAspectHeight 600
screen height 798
generate_aspect_width iAspectWidth 600
generate_aspect_height iAspectHeight 600
iAspectX 1 iAspectY 1
aspect width 798 aspect height 798
screen width 1440 map width 1026
sUserIP 128.111.61.94
sFingerprint a58dc743dabd82f3593ceabfefb3a2e9
fingerprint: a58dc743dabd82f3593ceabfefb3a2e9
ip: 128.111.61.94
userhostname: 128-111-61-94.vpn.ucsb.edu
protocol: http:
hostname: 127.0.0.1
pathname: /
port: 3989
search:
queries:

input$zoomtoextent 0
input$zoomtoprev 0
change database start Tasmania
sSelectDb Tasmania
sMarxanDir /Users/bbest/marxanui//marxan//Tasmania
ChangeDatabase start
sBLM 0.05
PrepareDisplay start
PrepareDisplay end
ChangeDatabase end
iAspectX 317746.78612611 iAspectY 342985.02329322
aspect width 739 aspect height 798
change database end Tasmania
m
iM 1
fEnableMap FALSE
fEnableLeaflet FALSE
blm
rblm 0.1
generate_aspect_width iAspectWidth 739
generate_aspect_height iAspectHeight 798
session$clientData
mrun
blm
rblm 0.05
generate_aspect_width iAspectWidth 739
generate_aspect_height iAspectHeight 798
```
