load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
	annfnames = (/"CNRM-CM5annual.txt",\
	"inmcm4annual.txt",\
	"MPI-ESM-LRannual.txt",\
	"NorESMannual.txt",\
	"mer_annual.txt",\
	"era_annual.txt",\
	"cfsr_annual.txt",\
	"narr_regrid_annual.txt"/)

	monfnames= (/"CNRM-CM5monthly.txt",\ 
	"inmcm4monthly.txt",\ 
	"MPI-ESM-LRmonthly.txt",\ 
	"NorESMmonthly.txt",\
	"mer_monthly.txt",\
	"era_monthly.txt",\
	"cfsr_monthly.txt",\
	"narr_regrid_monthly.txt"/)

	yrStart = 1979
	yrLast = 2005
	nyears = 31 ;some reanalysis uses 31 years, while most CMIP5 data have 26 years
	;nq = (yrLast-yrStart+1)*12



	nfiles = dimsizes(annfnames)
	annual = new((/nfiles,nyears/),float)
	annual@_FillValue = default_fillvalue("float")
	monthly = new((/nfiles,12/),float)

	do ifile = 0, nfiles-1
		ann = asciiread(annfnames(ifile),-1,"float")
		annual(ifile,0:dimsizes(ann)-1) = ann
		monthly(ifile,:) = asciiread(monfnames(ifile),-1,"float")
		print(ann)
		delete(ann)

	end do
	print(annual)
	;order is merra, era, cfsr, narr
	;		

;************************************************
; plotting parameters
;************************************************
 wks_type = "png"
 wks_type@wkWidth = 3200
 wks_type@wkHeight = 3200
 wks   = gsn_open_wks (wks_type,"mfluxwith1000")                ; open workstation
 gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")

 plot = new(1,graphic)

	xyLineColors = (/	27,\
						59,\
						11,171,\
						2,2,2,2/)
	xyMarkerColors = (/27,\
						59,\
						43,171,\
						18,82,178,147/)
	xyLineThicknesses = (/2,\
							2,\
							2,2,\
							3.5,3.5,3.5,3.5/)
	xyLineThicknesses = xyLineThicknesses*3
	xyMarkers = (/1,\
					1,\
					6,6,\
					16,16,16,16/)
	xyDashPatterns = (/	1,\
						1,\
						1,1,\
						1,2,3,0/)
	xyExplicitLegendLabels = (/	"CNRM-CM5",\
								"INMCM4", \
								"MPI","NorESM",\
								"MERRA","ERA","CFSR","NARR"/)


 res                  = True                     ; plot mods desired
 res@tiYAxisString   = "Moisture Flux (mm/day)"
 res@tiXAxisString   = ""
 res@gsnDraw              = True           ; don't draw
 res@gsnFrame             = True           ; don't advance frame
 res@gsnMaximize = True
 res@tmXBMode             = "Explicit"      ; label independently
 res@tmXBValues           = ispan(0,11,1)
 res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S",\
  "O","N","D"/)
  ;res@xyLineThicknessF = 3.0       ; double the thickness of the line

 time = ispan(0,11,1)
 
 res@tiMainString     = "Monthly Climatology"

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
  itemsIndex = ispan(0,nfiles-1,1)
  itemsIndex = itemsIndex(::-1) ;reverse
  res@lgItemOrder = itemsIndex
 plot  = gsn_csm_xy (wks,time,monthly,res) ; create plot

 delete(time)
 time = ispan(1979,2009,1)
 res@tmXBMode             = "Automatic" 
 res@tiMainString     = "Annual Climatology"
 plot  = gsn_csm_xy (wks,time,annual,res) ; create plot
 exit


	;cmip5annual = annual(:10,:)
	;reanlannual = annual(11:,:)
	;cmip5month = annual(:10,:)
	;reanlmonth = annual(11:,:)
	;5 here is 4 reanalysis + 1 CMIP5 on the same plot
	plotannual = new((/5,nyears/), float)
	plotmonthly = new((/5,12/), float)

	xyLineColors5 = new(5, integer)
	xyMarkerColors5 = new(5, integer)
	xyLineThicknesses5 = new(5, float)
	xyMarkers5 = new(5, integer)
	xyDashPatterns5 = new(5, integer)
	xyExplicitLegendLabels5 = new(5, string)
	

	delete(res@lgItemOrder)
	nAR5models = 12
	do i=0,11
		delete(res@xyLineColors)
		delete(res@xyMarkerColors)
		delete(res@xyLineThicknesses)
		delete(res@xyMarkers)
		delete(res@xyDashPatterns)
		delete(res@xyExplicitLegendLabels)
		delete(itemsIndex)
		delete(time)
		plotannual(0,:) = annual(i,:)
		plotannual(1:4,:) = annual(nAR5models:,:)
		plotmonthly(0,:) = monthly(i,:)
		plotmonthly(1:4,:) = monthly(nAR5models:,:)

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

		res@tiMainString = "Monthly Climatology"
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

		delete(time)
		res@tmXBMode             = "Automatic" 
		res@tiMainString     = "Annual Climatology"
		time = ispan(1979,2009,1)
		plot  = gsn_csm_xy (wks,time,plotannual,res) ; create plot


		;delete(xyLineColors)
		;delete(xyMarkerColors5)
		;delete(xyLineThicknesses5)
		;delete(xyMarkers5)
		;delete(xyDashPatterns5)
		;delete(xyExplicitLegendLabels5)
	end do


 







end
