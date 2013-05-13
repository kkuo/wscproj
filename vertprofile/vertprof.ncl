
; save new output variable; Author: Kevin Kuo, Allison Steiner
; ALS modifications 12 Oct 12 based on Zangvil et al., 2004


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------
; This function takes a regular XY plot with a linear Y axis, and 
; given a new set of Y axis values, changes the scale of the Y 
; axis to an irregular scale. 
;----------------------------------------------------------------------
undef("make_yaxis_irregular")
function make_yaxis_irregular(wks,plot,yvals)
local xaxis_type, xmin, xmax, ymin, ymax, yreverse
begin
;
; First retrieve some resource values from original plot so we
; can make sure new overlay object has same values.
;
; Depending on what axis you want to make "irregular", you may
; need to change which resources are retrieved and/or set.
;
  getvalues plot
    "trXAxisType"    : xaxis_type
    "trXMinF"        : xmin
    "trXMaxF"        : xmax
    "trYMinF"        : ymin
    "trYMaxF"        : ymax
    "trXReverse"     : xreverse
    "trYReverse"     : yreverse
  end getvalues

;
; Create a blank plot. This will enable us to make the Y axis
; of our existing plot irregular by overlaying our plot on the
; irregular object.
;
  bres                = True
  bres@gsnDraw        = False
  bres@trXAxisType    = xaxis_type
  bres@trXMinF        = xmin
  bres@trXMaxF        = xmax
  bres@trYMinF        = ymin
  bres@trYMaxF        = ymax
  bres@trXReverse     = xreverse
  bres@trYReverse     = yreverse

;---These two resources will change the Y axis to an irregular one.
  bres@trYAxisType    = "IrregularAxis"
  bres@trYCoordPoints = yvals

  blank_plot = gsn_csm_blank_plot(wks,bres)

;
; Overlay original plot on irregular object. This will cause a
; transformation to take place, and Y axis will be irregular.
;
  overlay(blank_plot,plot)

    return(blank_plot)
end

begin
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 0: Set constants and variables
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	g = 9.8        ; gravitational constant (m/s2)
	rhowat = 1000. ; density of water (kg/m3)

	; define GL region
	; need to be in float for gc_latlon
	slat = 40.
	nlat = 50.
	elon = 285.
	wlon = 265.

	;slat = 41.
	;nlat = 49.
	;elon = 285.;-360
	;wlon = 266.;-360
	;nmonths = 28*12 ; 1979 to 2007

	yrStart = 1979
	yrLast = 2005
	nq = (yrLast-yrStart+1)*12

	asciidir = "."
	asciifile = "VertProfile"
	asciiext = ".txt"

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 1:  Read in data and get dimensions
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

;	asciifnames = (/"/share/kkuo/AR5/scripts/vertprofile/CCSM4",\
;					"/share/kkuo/AR5/scripts/vertprofile/CNRM-CM5",\
;					"/share/kkuo/AR5/scripts/vertprofile/GFDL-ESM2G", \
;					"/share/kkuo/AR5/scripts/vertprofile/inmcm4",\
;					"/share/kkuo/AR5/scripts/vertprofile/IPSL-CM5A-LR",\
;					"/share/kkuo/AR5/scripts/vertprofile/HadESM2ES",\
;					"/share/kkuo/AR5/scripts/vertprofile/MIROC5",\
;					"/share/kkuo/AR5/scripts/vertprofile/MIROC-ESM",\
;					"/share/kkuo/AR5/scripts/vertprofile/MIROC-ESM-CHEM",\
;					"/share/kkuo/AR5/scripts/vertprofile/MPI-ESM-LR",\
;					"/share/kkuo/AR5/scripts/vertprofile/MRI-CGCM3",\
;					"/share/kkuo/AR5/scripts/vertprofile/NorESM"/)
;
	husfnames = (/	"/share/kkuo/AR5/CCSM4/keep/hus_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/hus_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/hus.nc", \
					"/share/kkuo/AR5/INMCM4/hus_Amon_inmcm4_1970_2005.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/hus_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/hus_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/hus_Amon_MIROC5_1970_2012.nc",\
					"/share/kkuo/AR5/MIROC-ESM/hus_Amon_MIROC-ESM_1850_2005.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/hus_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/hus_1970_to_2005.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/hus_Amon_197001-200512.nc",\
					"/share/kkuo/AR5/NorESM/hus_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

	uafnames = (/	"/share/kkuo/AR5/CCSM4/keep/ua_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/ua_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/ua.nc", \
					"/share/kkuo/AR5/INMCM4/ua_Amon_inmcm4_1970_2005.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/ua_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/ua_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/ua_Amon_MIROC5_1970_2012.nc",\
					"/share/kkuo/AR5/MIROC-ESM/ua_Amon_MIROC-ESM_1850_2005.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/ua_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/ua_1970_to_2005.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/ua_Amon_197001-200512.nc",\
					"/share/kkuo/AR5/NorESM/ua_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

	vafnames = (/	"/share/kkuo/AR5/CCSM4/keep/va_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/va_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/va.nc", \
					"/share/kkuo/AR5/INMCM4/va_Amon_inmcm4_1970_2005.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/va_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/va_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/va_Amon_MIROC5_1970_2012.nc",\
					"/share/kkuo/AR5/MIROC-ESM/va_Amon_MIROC-ESM_1850_2005.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/va_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/va_1970_to_2005.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/va_Amon_197001-200512.nc",\
					"/share/kkuo/AR5/NorESM/va_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

;	husfnames = (/	"/Volumes/Shared/AR5/GFDL-ESM2G/hus.nc", \ 
;					"/Volumes/Shared/AR5/HadESM2ES/hus_195912-200511.nc",\ 
;					"/Volumes/Shared/AR5/MPI-ESM-LR/hus_1970_to_2005.nc"/)
;
;	uafnames = (/	"/Volumes/Shared/AR5/GFDL-ESM2G/ua.nc", \ 
;					"/Volumes/Shared/AR5/HadESM2ES/ua_195912-200511.nc",\ 
;					"/Volumes/Shared/AR5/MPI-ESM-LR/ua_1970_to_2005.nc"/)
;
;	vafnames = (/	"/Volumes/Shared/AR5/GFDL-ESM2G/va.nc", \ 
;					"/Volumes/Shared/AR5/HadESM2ES/va_195912-200511.nc",\ 
;					"/Volumes/Shared/AR5/MPI-ESM-LR/va_1970_to_2005.nc"/)

	nfiles = dimsizes(husfnames)

	rad = 4.0*atan(1.0)/180.0
	re = 6371220.0
	rr = re*rad

	AllModelsW = new((/nfiles,25/),float)
	AllModelsE = new((/nfiles,25/),float)
	AllModelsN = new((/nfiles,25/),float)
	AllModelsS = new((/nfiles,25/),float)
	plev25 = (/10,20,30,100,200,300,500,700,1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,60000,70000,85000,92500,100000/)


	do ifile = 0, nfiles-1
		;if(isfilepresent(monFileStrings) .and. isfilepresent(annFileStrings))
		;	continue
		;end if

		;print(husfnames(ifile))
		;cmdstring = (/"ncdump -h ",husfnames(ifile)," | grep units"/)
		;cmd = str_concat(cmdstring)
		;system(cmd)
		;;print(husfnames(ifile))
		;cmdstring = (/"ncdump -h ",uafnames(ifile)," | grep units"/)
		;cmd = str_concat(cmdstring)
		;system(cmd)
		;;print(husfnames(ifile))
		;cmdstring = (/"ncdump -h ",vafnames(ifile)," | grep units"/)
		;cmd = str_concat(cmdstring)
		;system(cmd)

		fq = addfile(husfnames(ifile),"r")
		fu = addfile(uafnames(ifile),"r")
		fv = addfile(vafnames(ifile),"r")

		YYYYMM = cd_calendar(fq->time,-1)
		iStart = ind(YYYYMM.eq.(yrStart*100+1))
		iLast = ind(YYYYMM.eq.(yrLast*100+12))
		;HadESM2-ES data is numbered differently. December appears to be missing from 2005, but it is really not
		if(ismissing(iLast))
			iLast = dimsizes(YYYYMM)-1
		end if

		nq = iLast-iStart+1
		spechum1 = fq->hus(time|iStart:iLast, plev|:, lat|:, lon|:)
		spechum = spechum1({time|:},{plev|:},{lat|slat-5:nlat+5},{lon|wlon-10:elon+10})
		delete(spechum1)	
		;printMinMax(spechum,True)

		uwind1 = fu->ua(time|iStart:iLast, plev|:, lat|:,lon|:)
		uwind = uwind1({time|:}, {plev|:}, {lat|slat-5:nlat+5},{lon|wlon-10:elon+10})
		delete(uwind1)
		;printMinMax(uwind,True)

		vwind1 = fv->va(time|iStart:iLast,plev|:,lat|:,lon|:)
		vwind = vwind1({time|:},{plev|:}, {lat|slat-5:nlat+5},{lon|wlon-10:elon+10})
		delete(vwind1)
		;printMinMax(vwind,True)

		dims = dimsizes(spechum)
		levs = dims(1)
		lat2 = dims(2)
		lon2 = dims(3)
		;printVarSummary(spechum)
		lats = spechum&lat
		lons = spechum&lon
		;print(levs)

		; save new output variable
		flux  = new((/nq/),float)
		year  = new((/nq/),integer)
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; Step 3: calc dx and dy for each grid cell 
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; define a 2D dimension for each grid box
		dx = new((/lat2,lon2/),float)
		dy = new((/lat2,lon2/),float)

		dy2 = new(lat2, float)
		dy2(1:lat2-2) = dble2flt(abs(lats(2:lat2-1)-lats(1:lat2-2))*rr)
		dy2(0) = dy2(1)              ;set up dummy values for the first/last row
		dy2(lat2-1) = dy2(lat2-2)

		; calculate dx
		templon = spechum&lon
		;print(templon)
		dlontmp = templon(1)-templon(0)
		;print(dlontmp)
		dlon = dlontmp*rr
		; dx is a function of latitude
		do j=0,lon2-1
			do i=1,lat2-2
				dx(i,j) = doubletofloat(dlon*cos(lats(i)*rad))
			end do
		end do
		dx(0,:) = dx(1,:)           ; set up dummy vars for first and last row
		dx(lat2-1,:) = dx(lat2-2,:)

		dx!0 = "lat"
		dx!1 = "lon"
		dx&lat = lats
		dx&lon = lons
		dy!0 = "lat"
		dy!1 = "lon"
		dy&lat = lats
		dy&lon = lons
		; conform dy to the same grid (dy does not change with dx)
		dy = conform_dims(dimsizes(dx),dy2,0)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; Step 4: calculate dp for each pressure level
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		;fl = addfile("shum.mon.mean.nc","r")
		;pres = fl->level
		pres       = spechum&plev
		plev       = pres
		plev@units = "Pa"
		ptop       = min(plev)
		psfc       = 1013.*100.           ; use dummy value until get the 2D fields
		psfc@units = "Pa"
		dp         = dpres_plevel(plev,psfc,ptop,0)      ; calculate pressure thicknesses
		dp_3d      = conform_dims(dimsizes(uwind),dp,1)  ; need to make the dp array into 3D (use conform)
		print(plev)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; Step 5:  Moisture Flux Divergence Calculation
			;   Use Green's theorem as in Zangvil et al., 2004, Eqn 2
		;   Inner integral(dl) - q*v*dl [kg H2O/kg air * m/s * m = kg H2O m2/(kg air s)] 
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

		; conform dims for the dx and dy to make it 3D 
		newdx = conform(spechum,dx,(/2,3/))
		newdy = conform(spechum,dy,(/2,3/))
		copy_VarMeta(spechum,newdx)
		copy_VarMeta(spechum,newdy)
		;printVarSummary(newdx)

		;print(newdx(0,{slat},{wlon:elon}))

		; copy metadata
		quflx = spechum								
		qvflx = spechum 
		; uwind moisture flux divergence (kg m2)/(kg s) at each pressure level
		; for uwind, use the dy length scale
		quflx = spechum*uwind*newdy				
		; vwind moisture flux divergence (kg m2)/(kg s) at each pressure level
		; for vwind, use the dx length scale
		qvflx = spechum*vwind*newdx				

		; NOTE: We use dim_sum in the mflux script because we are summin
		; over the integrated result. Here, we should change to average 
		; because we are not performing an integration.
        WTimeLev = dim_avg_Wrap(quflx(:,:,{slat:nlat},{wlon}))
        ETimeLev = dim_avg_Wrap(quflx(:,:,{slat:nlat},{elon}))
        NTimeLev = dim_avg_Wrap(qvflx(:,:,{nlat},{wlon:elon}))
        STimeLev = dim_avg_Wrap(qvflx(:,:,{slat},{wlon:elon}))

		W = dim_avg_n_Wrap(WTimeLev,0)
		E = dim_avg_n_Wrap(ETimeLev,0)
		N = dim_avg_n_Wrap(NTimeLev,0)
		S = dim_avg_n_Wrap(STimeLev,0)

		print(W)
		do i = 0,25-1
			do j=0,dimsizes(plev)-1
				if (plev25(i) .eq. plev(j))
					AllModelsW(ifile,i) = W(j)
					AllModelsE(ifile,i) = E(j)
					AllModelsN(ifile,i) = N(j)
					AllModelsS(ifile,i) = S(j)
				end if
			end do
		end do
		print(AllModelsW)
		


		delete(YYYYMM)
		delete(iStart)
		delete(iLast)
		delete(spechum)
		delete(uwind)
		delete(vwind)
		delete(dims)
		delete(levs)
		delete(lat2)
		delete(lon2)
		delete(lats)
		delete(lons)
		delete(flux)
		delete(year)
		delete(dx)
		delete(dy)
		delete(dy2)
		delete(templon)
		delete(pres)
		delete(plev)
		delete(ptop)
		delete(psfc)
		delete(dp)
		delete(dp_3d)
		delete(newdx)
		delete(newdy)
		delete(quflx)
		delete(qvflx)
		delete(WTimeLev)
		delete(ETimeLev)
		delete(NTimeLev)
		delete(STimeLev)
		delete(W)
		delete(E)
		delete(N)
		delete(S)

	end do


;************************************************
; plotting parameters
;************************************************
	wks_type = "pdf"
	;wks_type@wkWidth=3200
	;wks_type@wkHeight=3200
	wks   = gsn_open_wks (wks_type,asciifile)                ; open workstation
	gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")

	plot = new(1,graphic)

	xyLineColorsList = (/11,27,43,\
						59,75,91,123,\
						123,123,11,155,171,\
						2,2,2,2/)
	xyMarkerColorsList = (/11,27,43,\
							59,75,91,20,\
							164,69,43,155,171,\
							18,82,178,147/)
	xyLineThicknessesList = (/2,2,2,\
								2,2,2,2,\
								2,2,2,2,2,\
								3.5,3.5,3.5,3.5/)
	xyLineThicknessesList = xyLineThicknessesList*3
	xyMarkersList = (/6,1,6,\
						1,6,1,16,\
						16,16,6,1,6,\
						16,16,16,16/)
	xyDashPatternsList = (/	1,1,1,\
							1,1,1,1,\
							1,1,1,1,1,\
							1,2,3,0/)
	xyExplicitLegendLabelsList = (/	"CCSM4","CNRM-CM5","GFDL-ESM2G",\
									"INMCM4", "IPSL","HadESM2-ES","MIROC5",\
									"MIROC-ESM","MIROC-ESM-CHEM","MPI","MRI","NorESM",\
									"MERRA","ERA","CFSR","NARR"/)
	;modelidx = (/2,5,10/)
	modelidx = ispan(0,nfiles-1,1)

	xyLineColors = xyLineColorsList(modelidx)
	xyMarkerColors = xyMarkerColorsList(modelidx)
	xyLineThicknesses = xyLineThicknessesList(modelidx) 
	xyMarkers = xyMarkersList(modelidx)
	xyDashPatterns = xyDashPatternsList(modelidx)
	xyExplicitLegendLabels = xyExplicitLegendLabelsList(modelidx)



	res					= True                     ; plot mods desired
	res@trYReverse		= True
	res@tiYAxisString   = "Pressure (hPa)"
	res@tiXAxisString   = "moisture flux (kg/kg*m*m/s)"
	res@gsnDraw			= False  
	res@gsnFrame		= False 
	res@gsnMaximize		= True
	;res@tmXBMode             = "Explicit"      ; label independently
	;res@tmXBValues           = ispan(0,11,1)
	;res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S",\
	;"O","N","D"/)
	res@xyLineThicknessF = 1.0       ; double the thickness of the line

	time = ispan(0,11,1)

	print(xyLineColors)

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
	res@pmLegendOrthogonalPosF = -0.005
	;res@pmLegendParallelPosF = 0.8
	res@lgPerimOn = True
	res@lgBoxMinorExtentF = 0.1
	res@pmLegendWidthF = 0.5
	res@pmLegendHeightF = 0.20
	itemsIndex = ispan(0,nfiles-1,1)
	itemsIndex = itemsIndex(::-1) ;reverse
	res@lgItemOrder = itemsIndex


	;print(W)
	;print(plev)

	irrypts = (/0,50,100,150,200,300,400,600,800,1000/)
	plev = plev25/100
	res@tiMainString = "West q*u"
	plot  = gsn_csm_xy (wks,AllModelsW,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	;maximize_output(wks,True)

	res@tiMainString = "North q*u"
	plot  = gsn_csm_xy (wks,AllModelsN,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	;maximize_output(wks,True)

	res@tiMainString = "East q*u"
	plot  = gsn_csm_xy (wks,AllModelsE,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	;maximize_output(wks,True)

	res@tiMainString = "South q*u"
	plot  = gsn_csm_xy (wks,AllModelsS,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	;maximize_output(wks,True)

	;writing out to ascii files and then reading them in does not make
	;sense because not all models treat levels the same way
	;asciiwrite(eastFileStrings,E)
	;asciiwrite(westFileStrings,W)
	;asciiwrite(southFileStrings,S)
	;asciiwrite(northFileStrings,N)


end
