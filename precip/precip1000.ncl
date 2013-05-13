load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
	g = 9.8        ; gravitational constant (m/s2)
	rhowat = 1000. ; density of water (kg/m3)

	; define GL region
	; need to be in float for gc_latlon
	slat = 40.
	nlat = 50.
	elon = 285.
	wlon = 265.

	;these are what Li uses
	;slat = 41.
	;nlat = 49.
	;elon = 285.;-360
	;wlon = 266.;-360

	yrStart = 1979
	yrLast = 2005
	nq = (yrLast-yrStart+1)*12

	asciidir = "."
	asciifile = "PrecipALL"
	asciiext = ".txt"


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 1:  Read in data and get dimensions
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
				;"/share/kkuo/AR5/BCC/pr_Amon_bcc-csm1-1_historical_r1i1p1_185001-201212.nc", \
				;"/share/kkuo/AR5/HadESM2AO/pr_Amon_HadGEM2-AO_historical_r1i1p1_186001-200512.nc",\
	fnames = (/	"/share/kkuo/AR5/CCSM4/keep/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/GFDL-ESM2G/pr_Amon_GFDL-ESM2G_historical_r1i1p1_194601-200512.nc", \ 
				"/share/kkuo/AR5/IPSL-CM5A-LR/pr_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/HadESM2ES/pr_195912-200511.nc",\
				"/share/kkuo/AR5/MIROC5/pr_Amon_MIROC5_historical_r1i1p1_185001-201212.nc",\
				"/share/kkuo/AR5/MIROC-ESM/pr_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/MIROC-ESM-CHEM/pr_Amon_MIROC-ESM-CHEM_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/MRI-CGCM3/pr_Amon_MRI-CGCM3_historical_r1i1p1_185001-200512.nc"/)

	nfiles = dimsizes(fnames)
	precip = new((/nfiles+4,12/),float)
	do ifile = 0, nfiles-1
		print(fnames(ifile))
		cmdstring = (/"ncdump -h ",fnames(ifile)," | grep pr:units"/)
		cmd = str_concat(cmdstring)
		system(cmd)


		fp = addfile(fnames(ifile),"r")
		YYYYMM = cd_calendar(fp->time,-1)
		iStart = ind(YYYYMM.eq.(yrStart*100+1))
		iLast = ind(YYYYMM.eq.(yrLast*100+12))
		if(ismissing(iLast))
			iLast = dimsizes(YYYYMM)-1
		end if

		pr1 = fp->pr(time|iStart:iLast, lat|:, lon|:)
		pr = pr1({time|:},{lat|slat:nlat},{lon|wlon:elon})
		dims_pr = dimsizes(pr)
		nmonths = dims_pr(0)
		prTimeAvg = dim_avg_n_Wrap(pr,(/1,2/))

		do imon=0,11
			precip(ifile,imon) = dim_avg(prTimeAvg(imon::12))
		end do
		delete(YYYYMM)
		delete(iStart)
		delete(iLast)
		delete(pr1)
		delete(pr)
		delete(prTimeAvg)
	end do
	precip = precip*86400 ; convert from kg/m^2/s to mm/day
	;order is merra, era, cfsr, narr
	;		
	precipObs = asciiread("/share/kkuo/AR5/scripts/precip/precip.txt",(/4,12/),"float")
	precip(nfiles:,:)=precipObs
	print(precip)

;************************************************
; plotting parameters
;************************************************
	wks_type = "png"
	wks_type@wkWidth=3200
	wks_type@wkHeight=3200
	wks   = gsn_open_wks (wks_type,"PrecipMissing1000")                ; open workstation
	gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")

	plot = new(1,graphic)

	xyLineColors = (/	11,43,\
						75,91,123,\
						123,123,155,\
						2,2,2,2/)
	xyMarkerColors = (/11,43,\
						75,91,20,\
						164,69,155,\
						18,82,178,147/)
	xyLineThicknesses = (/2,2,\
							2,2,2,\
							2,2,2,\
							3.5,3.5,3.5,3.5/)
	xyLineThicknesses = xyLineThicknesses*3
	xyMarkers = (/6,6,\
					6,1,16,\
					16,16,1,\
					16,16,16,16/)
	xyDashPatterns = (/	1,1,\
						1,1,1,\
						1,1,1,\
						1,2,3,0/)
	xyExplicitLegendLabels = (/	"CCSM4","GFDL-ESM2G",\
								"IPSL","HadESM2-ES","MIROC5",\
								"MIROC-ESM","MIROC-ESM-CHEM","MRI",\
								"MERRA","ERA","CFSR","NARR"/)



 res                  = True                     ; plot mods desired
 res@tiYAxisString   = "Precipitation (mm/day)"
 res@tiXAxisString   = ""
 res@gsnDraw              = True           ; don't draw
 res@gsnFrame             = True           ; don't advance frame
 res@gsnMaximize = True
 res@tmXBMode             = "Explicit"      ; label independently
 res@tmXBValues           = ispan(0,11,1)
 res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S",\
  "O","N","D"/)
  res@xyLineThicknessF = 3.0       ; double the thickness of the line

 time = ispan(0,11,1)
 
 res@tiMainString     = asciifile+" lat=40N-50N,lon=95W-75W"

  res@xyLineColors = xyLineColors
  res@xyMarkerColors = xyMarkerColors
  res@xyLineThicknesses = xyLineThicknesses
  res@xyMarkLineMode = "MarkLines"
  res@xyMarkerSizeF = 0.008
  res@xyMarkers = xyMarkers
  res@xyDashPatterns = xyDashPatterns
  res@pmLegendDisplayMode = "Always"
  res@xyExplicitLegendLabels = xyExplicitLegendLabels
  ;res@trYMinF = -20
  ;res@pmLegendWidthF = 0.12
  ;res@pmLegendHeightF = 0.15
  res@pmLegendSide = "Bottom"
  res@pmLegendOrthogonalPosF = -0.07
  ;res@pmLegendParallelPosF = 0.8
  res@lgPerimOn = False
  res@lgBoxMinorExtentF = 0.3
  res@pmLegendWidthF = 0.5
  ;res@pmLegendHeightF = 1.5
  itemsIndex = ispan(0,nfiles+4-1,1)
  itemsIndex = itemsIndex(::-1) ;reverse
  res@lgItemOrder = itemsIndex
 plot  = gsn_csm_xy (wks,time,precip,res) ; create plot

	plotmonthly = new((/5,12/), float)

	xyLineColors5 = new(5, integer)
	xyMarkerColors5 = new(5, integer)
	xyLineThicknesses5 = new(5, float)
	xyMarkers5 = new(5, integer)
	xyDashPatterns5 = new(5, integer)
	xyExplicitLegendLabels5 = new(5, string)
	

	delete(res@lgItemOrder)
	nAR5models = 8
	do i=0,nAR5models-1
		delete(res@xyLineColors)
		delete(res@xyMarkerColors)
		delete(res@xyLineThicknesses)
		delete(res@xyMarkers)
		delete(res@xyDashPatterns)
		delete(res@xyExplicitLegendLabels)
		delete(itemsIndex)
		plotmonthly(0,:) = precip(i,:)
		plotmonthly(1:4,:) = precip(nAR5models:,:)

		xyLineColors5(0) = xyLineColors(i)
		xyLineColors5(1:4) = xyLineColors(nAR5models:)
		xyMarkerColors5(0) = xyMarkerColors(i)
		xyMarkerColors5(1:4) = xyMarkerColors(nAR5models:)
		xyLineThicknesses5(0) = xyLineThicknesses(i)
		xyLineThicknesses5(1:4) = xyLineThicknesses(nAR5models:)
		xyMarkers5(0) = xyMarkers(i)
		xyMarkers5(1:4) = xyMarkers(nAR5models:)
		xyDashPatterns5(0) = xyDashPatterns(i)
		xyDashPatterns5(1:4) = xyDashPatterns(nAR5models:)
		xyExplicitLegendLabels5(0) = xyExplicitLegendLabels(i)
		xyExplicitLegendLabels5(1:4) = xyExplicitLegendLabels(nAR5models:)

		res@tiMainString = "Precip Monthly " + xyExplicitLegendLabels(i)
		res@tmXBMode = "Explicit"      ; label independently
		res@tmXBValues = ispan(0,11,1)
		res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
		res@xyLineColors = xyLineColors5
		res@xyMarkerColors = xyMarkerColors5
		res@xyLineThicknesses = xyLineThicknesses5
		res@xyMarkLineMode = "MarkLines"
		res@xyMarkerSizeF = 0.008
		res@xyMarkers = xyMarkers5
		res@xyDashPatterns = xyDashPatterns5
		res@pmLegendDisplayMode = "Always"
		res@xyExplicitLegendLabels = xyExplicitLegendLabels5


		res@pmLegendSide = "Bottom"
		res@pmLegendOrthogonalPosF = -0.07
		;res@pmLegendParallelPosF = 0.8
		res@lgPerimOn = False
		res@lgBoxMinorExtentF = 0.3
		res@pmLegendWidthF = 0.5
		;res@pmLegendHeightF = 1.5
		itemsIndex = ispan(0,4,1)
		itemsIndex = itemsIndex(::-1) ;reverse
		res@lgItemOrder = itemsIndex
		time = ispan(0,11,1)
		plot  = gsn_csm_xy (wks,time,plotmonthly,res) ; create plot


		;delete(xyLineColors)
		;delete(xyMarkerColors5)
		;delete(xyLineThicknesses5)
		;delete(xyMarkers5)
		;delete(xyDashPatterns5)
		;delete(xyExplicitLegendLabels5)
	end do



end
