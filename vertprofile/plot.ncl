; save new output variable; Author: Kevin Kuo, Allison Steiner
; ALS modifications 12 Oct 12 based on Zangvil et al., 2004

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../common/modelnames.ncl"
load "make_yaxis_irregular.ncl"


begin

	if(.not. isvar("bExecuteOnEdison")) then
		bExecuteOnEdison = True
		sPathPrefix= "/share/kkuo/AR5/"
	else
		sPathPrefix= "/Volumes/Shared/AR5/kkuo/AR5/"
	end if

	if(.not. isvar("nModels")) then
		nModels = getmodelcount()
	end if

	AllModelsW = new((/nModels,25/),float)
	AllModelsE = new((/nModels,25/),float)
	AllModelsN = new((/nModels,25/),float)
	AllModelsS = new((/nModels,25/),float)
	plev25 = (/10,20,30,100,200,300,500,700,1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,60000,70000,85000,92500,100000/)

	do idx = 0,nModels-1
		modelNameE = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(idx) + "_E.txt"
		modelNameW = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(idx) + "_W.txt"
		modelNameS = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(idx) + "_S.txt"
		modelNameN = sPathPrefix + "scripts/vertprofile/" + getmodelnamebyindex(idx) + "_N.txt"

		E1 = asciiread(modelNameE,-1,"float")
		W1 = asciiread(modelNameW,-1,"float")
		S1 = asciiread(modelNameS,-1,"float")
		N1 = asciiread(modelNameN,-1,"float")

		length = dimsizes(E1)/2
		E2 = onedtond(E1,(/length,2/))
		W2 = onedtond(W1,(/length,2/))
		S2 = onedtond(S1,(/length,2/))
		N2 = onedtond(N1,(/length,2/))

		plev = E2(:,0)
		E = E2(:,1)
		W = W2(:,1)
		S = S2(:,1)
		N = N2(:,1)
		E@_FillValue = -9999
		W@_FillValue = -9999
		S@_FillValue = -9999
		N@_FillValue = -9999

		print(plev)
		print(E)
		
		
		; basically, some models have the full 25 levels while others have a subset
		; this is to snure those models that do not have the full 25 levels
		; can still be plotted by matching the correct plevs from the subset
		; to the full set and let the unmatched plev[i] be missing
		do i = 0,25-1
			do j=0,dimsizes(plev)-1
				if (plev25(i) .eq. plev(j))
					AllModelsW(idx,i) = W(j)
					AllModelsE(idx,i) = E(j)
					AllModelsN(idx,i) = N(j)
					AllModelsS(idx,i) = S(j)
				end if
			end do
		end do
		delete(E1)
		delete(W1)
		delete(S1)
		delete(N1)
		delete(E2)
		delete(W2)
		delete(S2)
		delete(N2)
		delete(E)
		delete(W)
		delete(S)
		delete(N)
		delete(plev)

	end do


		
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
	xyLineThicknessesList = xyLineThicknessesList
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
	;modelidx = (/1,5,10/)
	modelidx = ispan(0,nModels-1,1)
	print(modelidx)
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
	;res@gsnMaximize		= True
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
	res@lgJustification = "TopRight"
	res@pmLegendOrthogonalPosF = -1.18
	res@pmLegendParallelPosF = 1.0
	res@lgPerimOn = True
	res@lgAutoManage = True
	;res@lgBoxMinorExtentF = 0.1
	;res@lgBoxMajorExtentF = 0.0
	res@lgLabelFontHeightF = 0.005
	res@pmLegendWidthF = 0.2
	res@pmLegendHeightF = 0.1
	itemsIndex = ispan(0,nModels-1,1)
	itemsIndex = itemsIndex(::-1) ;reverse
	res@lgItemOrder = itemsIndex

	;res@gsnPaperMargin = 0.3
	;res@gsnPanelBottom = 0.4


	;print(W)
	;print(plev)

	irrypts = (/0,50,100,150,200,300,400,600,800,1000/)
	plev = plev25/100

	res@tiMainString = "East q*u"
	plot  = gsn_csm_xy (wks,AllModelsE,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,True)

	res@tiMainString = "West q*u"
	plot  = gsn_csm_xy (wks,AllModelsW,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,True)

	res@tiMainString = "South q*u"
	plot  = gsn_csm_xy (wks,AllModelsS,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,True)

	res@lgJustification = "TopLeft"
	res@pmLegendOrthogonalPosF = -1.18
	res@pmLegendParallelPosF = 0.0

	res@tiMainString = "North q*u"
	plot  = gsn_csm_xy (wks,AllModelsN,plev,res) ; create plot
	irregular_plot = make_yaxis_irregular(wks,plot,irrypts)
	maximize_output(wks,True)

	;writing out to ascii files and then reading them in does not make
	;sense because not all models treat levels the same way
	;asciiwrite(eastFileStrings,E)
	;asciiwrite(westFileStrings,W)
	;asciiwrite(southFileStrings,S)
	;asciiwrite(northFileStrings,N)


end
