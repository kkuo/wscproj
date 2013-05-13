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
	asciifile = "P minus E"
	asciiext = ".txt"


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 1:  Read in data and get dimensions
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

				;"/share/kkuo/AR5/BCC/evap_OImon_bcc-csm1-1_historical_r1i1p1_185001-201212.nc" 
				;"/share/kkuo/AR5/HadESM2AO/pr_Amon_HadGEM2-AO_historical_r1i1p1_186001-200512.nc",\
	fnamesPR = (/	"/share/kkuo/AR5/CCSM4/keep/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/pr_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/pr_Amon_GFDL-ESM2G_historical_r1i1p1_194601-200512.nc", \ 
					"/share/kkuo/AR5/INMCM4/pr_Amon_inmcm4_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/pr_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/pr_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/pr_Amon_MIROC5_historical_r1i1p1_185001-201212.nc",\
					"/share/kkuo/AR5/MIROC-ESM/pr_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/pr_Amon_MIROC-ESM-CHEM_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/pr_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/pr_Amon_MRI-CGCM3_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/NorESM/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc"/)

	fnamesEV = (/	"/share/kkuo/AR5/CCSM4/keep/evspsbl_Amon_CCSM4_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/evspsbl_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/evspsbl_196101-200512.nc", \
					"/share/kkuo/AR5/INMCM4/evspsbl_Amon_inmcm4_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/evspsbl_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/evspsbl_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/evspsbl_Amon_MIROC5_historical_r1i1p1_185001-201212.nc",\
					"/share/kkuo/AR5/MIROC-ESM/evspsbl_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/evspsbl_Amon_MIROC-ESM-CHEM_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/evspsbl_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/evspsbl_Amon_MRI-CGCM3_historical_r1i1p1_185001-200512.nc",\
					"/share/kkuo/AR5/NorESM/evspsbl_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc"/)


	nAR5models = dimsizes(fnamesPR)
	nObs = 4
	nfiles = nAR5models + nObs
	evap = new((/nfiles,12/),float)
	precip = new((/nfiles,12/),float)
	;evap = new((/nfiles,12/),float)
	;precip = new((/nfiles,12/),float)
	do ifile = 0, nAR5models-1
		print(fnamesEV(ifile))
		cmdstring = (/"ncdump -h ",fnamesEV(ifile)," | grep pr:units"/)
		cmd = str_concat(cmdstring)
		system(cmd)

		fp = addfile(fnamesEV(ifile),"r")
		YYYYMM = cd_calendar(fp->time,-1)
		iStart = ind(YYYYMM.eq.(yrStart*100+1))
		iLast = ind(YYYYMM.eq.(yrLast*100+12))
		if(ismissing(iLast))
			iLast = dimsizes(YYYYMM)-1
		end if
		print(iStart)
		print(yrStart*100)

		ev1 = fp->evspsbl(time|iStart:iLast, lat|:, lon|:)
		ev = ev1({time|:},{lat|slat:nlat},{lon|wlon:elon})

		dims_ev = dimsizes(ev)
		nmonths = dims_ev(0)
		evTimeAvg = dim_avg_n_Wrap(ev,(/1,2/))

		do imon=0,11
			evap(ifile,imon) = dim_avg(evTimeAvg(imon::12))
		end do
		delete(YYYYMM)
		delete(iStart)
		delete(iLast)
		delete(ev1)
		delete(ev)
		delete(evTimeAvg)
		delete(fp)




		print(fnamesPR(ifile))
		cmdstring = (/"ncdump -h ",fnamesPR(ifile)," | grep pr:units"/)
		cmd = str_concat(cmdstring)
		system(cmd)
		

		fp = addfile(fnamesPR(ifile),"r")
		YYYYMM = cd_calendar(fp->time,-1)
		iStart = ind(YYYYMM.eq.(yrStart*100+1))
		iLast = ind(YYYYMM.eq.(yrLast*100+12))
		if(ismissing(iLast))
			iLast = dimsizes(YYYYMM)-1
		end if
		print(iStart)
		print(yrStart*100)

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
	evap = evap*86400 ; convert from kg/m^2/s to mm/day
	precip = precip*86400 ; convert from kg/m^2/s to mm/day

	;order is merra, era, cfsr, narr
	;		
	evapObs = asciiread("/share/kkuo/AR5/scripts/evap/evap.txt",(/nObs,12/),"float")
	precipObs = asciiread("/share/kkuo/AR5/scripts/precip/precip.txt",(/nObs,12/),"float")
	evap(nAR5models:,:)=evapObs
	precip(nAR5models:,:)=precipObs
	PmE = precip - evap
	;print(evap)
	;print(precip)

	; process mflux

	monfnames= (/"/share/kkuo/AR5/scripts/mflux/CCSM4monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/CNRM-CM5monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/GFDL-ESM2Gmonthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/HadESM2ESmonthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/inmcm4monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/IPSL-CM5A-LRmonthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/MIROC5monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/MIROC-ESMmonthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/MIROC-ESM-CHEMmonthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/MPI-ESM-LRmonthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/MRI-CGCM3monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/NorESMmonthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/mer_monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/era_monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/cfsr_monthly.txt",\
	"/share/kkuo/AR5/scripts/mflux/narr_regrid_monthly.txt"/)

	monthly = new((/nfiles,12/),float)

	do ifile = 0, nfiles-1
		monthly(ifile,:) = asciiread(monfnames(ifile),-1,"float")
	end do



;************************************************
; plotting parameters
;************************************************
	wks_type = "png"
	wks_type@wkWidth=3200
	wks_type@wkHeight=3200
	wks   = gsn_open_wks (wks_type,"mflux_PE_panel")                ; open workstation
	gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")

	plot = new(2,graphic)

	xyLineColors = (/	11,27,43,\
			59,75,91,123,\
			123,123,11,155,171,\
			2,2,2,2/)
	xyMarkerColors = (/11,27,43,\
			59,75,91,20,\
			164,69,43,155,171,\
			18,82,178,147/)
	xyLineThicknesses = (/2,2,2,\
			2,2,2,2,\
			2,2,2,2,2,\
			3.5,3.5,3.5,3.5/)
	xyLineThicknesses = xyLineThicknesses*3
	xyMarkers = (/6,1,6,\
			1,6,1,16,\
			16,16,6,1,6,\
			16,16,16,16/)
	xyDashPatterns = (/	1,1,1,\
			1,1,1,1,\
			1,1,1,1,1,\
			1,2,3,0/)
	xyExplicitLegendLabels = (/	"CCSM4","CNRM-CM5","GFDL-ESM2G",\
			"INMCM4", "IPSL","HadESM2-ES","MIROC5",\
			"MIROC-ESM","MIROC-ESM-CHEM","MPI","MRI","NorESM",\
			"MERRA","ERA","CFSR","NARR"/)




	res                  = True                     ; plot mods desired
	res@tiYAxisString   = "Difference (mm/day)"
	res@tiXAxisString   = ""
	res@gsnDraw              = False           ; don't draw
	res@gsnFrame             = False           ; don't advance frame
	res@gsnMaximize = False
	res@tmXBMode             = "Explicit"      ; label independently
	res@tmXBValues           = ispan(0,11,1)
	res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S",\
			"O","N","D"/)

	time = ispan(0,11,1)

	;res@tiMainString     = asciifile+" lat=40N-50N,lon=95W-75W"

	res@xyLineColors = xyLineColors
	res@xyMarkerColors = xyMarkerColors
	res@xyLineThicknesses = xyLineThicknesses
	res@xyMarkLineMode = "MarkLines"
	res@xyMarkerSizeF = 0.008
	res@xyMarkers = xyMarkers
	res@xyDashPatterns = xyDashPatterns
	res@pmLegendDisplayMode = "Never"
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
	itemsIndex = ispan(0,nfiles-1,1)
	itemsIndex = itemsIndex(::-1) ;reverse
	res@lgItemOrder = itemsIndex
	res@gsnLeftString = "Moisture Flux (mm/day)"
	plot(0)  = gsn_csm_xy (wks,time,monthly,res) ; create plot
	res@gsnLeftString = "P-E (mm/day)"
	plot(1)  = gsn_csm_xy (wks,time,PmE,res) ; create plot

	resP = True
	resP@gsnMaximize = True
	resP@txString = "Moisture Flux and P-E, All Models"
	gsn_panel(wks,plot,(/1,2/),resP)




	nLines = 1
	plotmonthly = new((/nLines,12,2/), float)

	xyLineColors5 = new(nLines, integer)
	xyMarkerColors5 = new(nLines, integer)
	xyLineThicknesses5 = new(nLines, float)
	xyMarkers5 = new(nLines, integer)
	xyDashPatterns5 = new(nLines, integer)
	xyExplicitLegendLabels5 = new(nLines, string)
	
	delete(res@lgItemOrder)

	do i=0,nfiles-1
		delete(res@xyLineColors)
		delete(res@xyMarkerColors)
		delete(res@xyLineThicknesses)
		delete(res@xyMarkers)
		delete(res@xyDashPatterns)
		delete(res@xyExplicitLegendLabels)
		;delete(itemsIndex)
		plotmonthly(0,:,0) = monthly(i,:)
		plotmonthly(0,:,1) = PmE(i,:)
		;plotmonthly(1:4,:) = PmE(nAR5models:,:)

		xyLineColors5(0) = xyLineColors(i)
		;xyLineColors5(1:4) = xyLineColors(nAR5models:)
		xyMarkerColors5(0) = xyMarkerColors(i)
		;xyMarkerColors5(1:4) = xyMarkerColors(nAR5models:)
		xyLineThicknesses5(0) = xyLineThicknesses(i)
		;xyLineThicknesses5(1:4) = xyLineThicknesses(nAR5models:)
		xyMarkers5(0) = xyMarkers(i)
		;xyMarkers5(1:4) = xyMarkers(nAR5models:)
		xyDashPatterns5(0) = xyDashPatterns(i)
		;xyDashPatterns5(1:4) = xyDashPatterns(nAR5models:)
		xyExplicitLegendLabels5(0) = xyExplicitLegendLabels(i)
		;xyExplicitLegendLabels5(1:4) = xyExplicitLegendLabels(nAR5models:)

		;res@tiMainString = "P-E Monthly " + xyExplicitLegendLabels(i) 
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
		res@pmLegendDisplayMode = "Never"
		res@xyExplicitLegendLabels = xyExplicitLegendLabels5


		res@pmLegendSide = "Bottom"
		res@pmLegendOrthogonalPosF = -0.05
		res@pmLegendParallelPosF = 0.5
		res@lgPerimOn = False
		res@lgBoxMinorExtentF = 0.3
		res@pmLegendWidthF = 0.2
		res@pmLegendHeightF = 0.066
		;itemsIndex = ispan(0,4,1)
		;itemsIndex = itemsIndex(::-1) ;reverse
		;res@lgItemOrder = itemsIndex
		time = ispan(0,11,1)
		res@gsnLeftString = "Moisture Flux (mm/day)"
		plot(0)  = gsn_csm_xy (wks,time,plotmonthly(:,:,0),res) ; create plot
		res@gsnLeftString = "P-E (mm/day)"
		plot(1)  = gsn_csm_xy (wks,time,plotmonthly(:,:,1),res) ; create plot
		resP@txString = xyExplicitLegendLabels(i) 
		gsn_panel(wks,plot,(/1,2/),False)


		;delete(xyLineColors)
		;delete(xyMarkerColors5)
		;delete(xyLineThicknesses5)
		;delete(xyMarkers5)
		;delete(xyDashPatterns5)
		;delete(xyExplicitLegendLabels5)
	end do




end
