; save new output variable; Author: Kevin Kuo, Allison Steiner
; ALS modifications 12 Oct 12 based on Zangvil et al., 2004

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../common/modelnames.ncl"

begin
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 0: Set constants and variables
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	g = 9.8        ; gravitational constant (m/s2)
	rhowat = 1000. ; density of water (kg/m3)

	; define GL region
	if(.not. isvar("bUseLiLatLon")) then
		bUseLiLatLon = True
	end if
	if (bUseLiLatLon .eq. True) then
		slat = 40.
		nlat = 50.
		elon = 285.
		wlon = 265.
	else
		slat = 41.
		nlat = 49.
		elon = 285.;-360
		wlon = 266.;-360
	end if

	yrStart = 1979
	yrLast = 2005
	nq = (yrLast-yrStart+1)*12


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 1:  Read in data and get dimensions
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	if(.not. isvar("bExecuteOnEdison")) then
		bExecuteOnEdison = True
		sPathPrefix= "/share/kkuo/AR5/"
	else
		sPathPrefix= "/Volumes/Shared/AR5/kkuo/AR5/"
	end if

	if(.not. isvar("iFileIndex")) then
		nfiles = getmodelcount()
		print("please pass in a file index when calling this file, for example, iFileIndex=N where N is less than "+nfiles)
		exit
	end if

	nfiles = dimsizes(husfnames)

	rad = 4.0*atan(1.0)/180.0
	re = 6371220.0
	rr = re*rad

;	AllModelsW = new((/nfiles,25/),float)
;	AllModelsE = new((/nfiles,25/),float)
;	AllModelsN = new((/nfiles,25/),float)
;	AllModelsS = new((/nfiles,25/),float)
;	plev25 = (/10,20,30,100,200,300,500,700,1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,60000,70000,85000,92500,100000/)

	husfname = sPathPrefix + gethusmodelbyindex(iFileIndex)
	vafname = sPathPrefix + getvamodelbyindex(iFileIndex)
	uafname = sPathPrefix + getuamodelbyindex(iFileIndex)
	fq = addfile(husfname,"r")
	fu = addfile(uafname,"r")
	fv = addfile(vafname,"r")

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
	quflx@_FillValue = -9999
	qvflx@_FillValue = -9999
	print(quflx@_FillValue)
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

	;print(WTimeLev)
	W = dim_avg_n_Wrap(WTimeLev,0)
	E = dim_avg_n_Wrap(ETimeLev,0)
	N = dim_avg_n_Wrap(NTimeLev,0)
	S = dim_avg_n_Wrap(STimeLev,0)

	sFileStringE = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(iFileIndex) + "_E.txt"
	sFileStringW = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(iFileIndex) + "_W.txt"
	sFileStringS = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(iFileIndex) + "_S.txt"
	sFileStringN = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(iFileIndex) + "_N.txt"

	print(sFileStringW)
	print(W)
	
	Eout = new ((/dimsizes(plev),2/), double)
	Wout = new ((/dimsizes(plev),2/), double)
	Nout = new ((/dimsizes(plev),2/), double)
	Sout = new ((/dimsizes(plev),2/), double)

	Eout(:,0) = plev
	Eout(:,1) = E
	Wout(:,0) = plev
	Wout(:,1) = W
	Nout(:,0) = plev
	Nout(:,1) = N 
	Sout(:,0) = plev
	Sout(:,1) = S 


	asciiwrite(sFileStringE,Eout)
	asciiwrite(sFileStringW,Wout)
	asciiwrite(sFileStringS,Sout)
	asciiwrite(sFileStringN,Nout)

	exit

		;print(W)
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
		;print(AllModelsW)
		


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



;************************************************
; plotting parameters
;************************************************
	wks_type = "pdf"
	;wks_type@wkWidth=3200
	;wks_type@wkHeight=3200
	wks   = gsn_open_wks (wks_type,"VertProfile")                ; open workstation
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
	modelidx = (/2,5,10/)
	;modelidx = ispan(0,nfiles-1,1)

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

	;print(xyLineColors)

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
	;res@pmLegendSide = "Bottom"
	res@lgJustification = "TopLeft"
	res@pmLegendOrthogonalPosF = -1.18
	res@pmLegendParallelPosF = 1.0
	res@lgPerimOn = True
	res@lgAutoManage = True
	;res@lgBoxMinorExtentF = 0.1
	;res@lgBoxMajorExtentF = 0.0
	res@lgLabelFontHeightF = 0.005
	res@pmLegendWidthF = 0.2
	res@pmLegendHeightF = 0.05
	itemsIndex = ispan(0,nfiles-1,1)
	itemsIndex = itemsIndex(::-1) ;reverse
	res@lgItemOrder = itemsIndex

	;res@gsnPaperMargin = 0.3
	;res@gsnPanelBottom = 0.4


	;print(W)
	;print(plev)

	irrypts = (/0,50,100,150,200,300,400,600,800,1000/)
	plev = plev25/100
	res@tiMainString = "West q*u"
	plot  = gsn_csm_xy (wks,AllModelsW,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,False)

	res@tiMainString = "North q*u"
	plot  = gsn_csm_xy (wks,AllModelsN,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,False)

	res@tiMainString = "East q*u"
	plot  = gsn_csm_xy (wks,AllModelsE,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,False)

	res@tiMainString = "South q*u"
	plot  = gsn_csm_xy (wks,AllModelsS,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,False)

	;writing out to ascii files and then reading them in does not make
	;sense because not all models treat levels the same way
	;asciiwrite(eastFileStrings,E)
	;asciiwrite(westFileStrings,W)
	;asciiwrite(southFileStrings,S)
	;asciiwrite(northFileStrings,N)


end
