clear
clc
close all
%% instructions
% The input file containing aspect data must be called 'aspect...' or another name that ensures it to be the first in alphabetical order. It must be in TIFF or GeoTIFF format.
% The inputfile containing slope data must be called 'slope..." or another name that ensures it to be the last in alphabetical order. It must be in TIFF or GeoTIFF format.
% The input file containing orientation of discontinuity sets must contain, in this exact order: [dip, dip direction, friction angle]. It must be in txt, csv or xlsx format.
% Files containing aspect and slope data must be placed in a different directory from the one containing orientation file.
% Output files are saved in DIR_output directory. They consist of TIFF or GeoTIFF files (same format as the one of the input slope and aspect files) that can be directly input in Gis environment. 

%% directories of input files
DIR_tiff = 'C:\Users\Gessica\Desktop\input\tiff\';% input data directory containing aspect and slope files 
DIR_orientation = 'C:\Users\Gessica\Desktop\input\csv\';% input data directory containing discontinuity sets orientation
DIR_output = 'C:\Users\Gessica\Desktop\output\';% output data directory, in which the results will be saved
filename = 'orientation.xlsx';% name of the file containing orientation data of the discontinuity sets
fileout2D = '2D_';% prefix of files containing the results of kinematic tests on planar (2D) sliding (it will be followed by the id of the considered set and 'sum' in the case of sum of all possible planar sliding kinematisms)
fileout3D = '3D_';% prefix of files containing the results of kinematic tests on wedge (3D) sliding (it will be followed by the id of the considered combination of sets and 'sum' in the case of sum of all possible wedge sliding kinematisms)
fileoutTOP = 'TOP_';% prefix of files containing the results of kinematic tests on flexural toppling (it will be followed by the id of the considered set and 'sum' in the case of sum of all possible flexural toppling kinematisms)
fileoutfinal = 'RESULTS';% name of the file containing the results of all the possible kinematisms
fileoutsource = 'SOURCES';% name of the file containing source areas

%% DO NOT CHANGE THE CODE FROM HERE ON OUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir(DIR_output);
IMlist = dir(DIR_tiff);
IMlist(1:2) = [];
IMlist(3) = [];
nImm = length(IMlist);
jImm = 1;

if nImm==2

    while (jImm <= nImm)
        if jImm == 1
            aspect = imread([DIR_tiff IMlist(jImm).name]);
            info_aspect = geotiffinfo([DIR_tiff IMlist(jImm).name]);
        else
            slope = imread([DIR_tiff IMlist(jImm).name]);
            info_slope = geotiffinfo([DIR_tiff IMlist(jImm).name]);
        end
        jImm = jImm+1;

    end

    %% orientation of discontinuity sets
    Orientazione = readtable([DIR_orientation filename]);
    Orientazione = table2array(Orientazione);
    nset = size(Orientazione,1);

    %% operations
    aspect(aspect==-9999)= NaN;
    slope(slope==-9999)= NaN;

    VER2Dsum = zeros(size(aspect,1),size(aspect,2));
    VER3Dsum = zeros(size(aspect,1),size(aspect,2));
    VERTOPsum = zeros(size(aspect,1),size(aspect,2));
    VERsum = zeros(size(aspect,1),size(aspect,2));

    %% MARKLAND'S TEST FOR PLANAR (2D) SLIDING

    for i = 1:nset

        % riclassification of aspect
        M_dipdir = Orientazione(i,2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - M_dipdir);
        Tmax = 40;
        Tmin = 0;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmax)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmin)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Tphi = Orientazione(i,3);
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<Tphi)= NaN;
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        Verifica2D = Condizione_orientazione_reclass & Condizione_pendenza;

        Verifica2D = double(Verifica2D);
        Verifica2Dnonriclass = Verifica2D;

        % write geotiff in output
        filenameout = [fileout2D num2str(i)];
        geotiffwrite([DIR_output filenameout],Verifica2D,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VER2D(:,:,i) = Verifica2D;
        VER2Dsum = VER2Dsum + Verifica2Dnonriclass;
    end

    filenameout = [fileout2D 'sum'];
    geotiffwrite([DIR_output filenameout],VER2Dsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);


    %% MARKLAND'S TEST FOR WEDGE (3D) SLIDING
    
    % identification of all the possible combinations of sets
    n_col = 2;
    t = 1:nset;
    y = repmat(t,1,nset^(n_col-1))';
    x = repmat(t,nset,nset^(n_col-2));
    xr = reshape(x,[size(y) 1]);
    C = [y xr];
    comb = sort(C,2);
    comb = unique(comb, 'rows');
    sy = [1:size(comb(:,1))]';
    ind = comb(:,1)-comb(:,2);
    ind(:,2) = sy;
    ind(ind==0,:) = [];
    Comb_set = comb(ind(:,2),:);

    for j = 1:size(Comb_set,1)
        dip1 = Orientazione(Comb_set(j,1),1);
        dipdir1 = Orientazione(Comb_set(j,1),2);
        dip2 = Orientazione(Comb_set(j,2),1);
        dipdir2 = Orientazione(Comb_set(j,2),2);
        n1 = [sin(dip1*pi/180)*sin(dipdir1*pi/180) sin(dip1*pi/180)*cos(dipdir1*pi/180) cos(dip1*pi/180)];
        n2 = [sin(dip2*pi/180)*sin(dipdir2*pi/180) sin(dip2*pi/180)*cos(dipdir2*pi/180) cos(dip2*pi/180)];
        l1 = norm(n1);
        l2 = norm(n2);
        pv = cross(n1/l1,n2/l2);

        if pv(3) < 0
            pv = -pv;
        end
        dip = acos(abs(pv(3))/norm(pv,'fro'))*180/pi;

        angle = acos(pv(2)/norm(pv(1:2),'fro'))*180/pi;
        if pv(3) >= 0
            if pv(1) >= 0
                dipDir = angle;
            else
                dipDir = 360-angle;
            end
        else
            if pv(1) >= 0
                dipDir = 180 + angle;
            else
                dipDir = 180 - angle;
            end
        end
        intersezione = [dip dipDir];

        % riclassification of aspect
        I_dipdir = intersezione(2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = aspect - I_dipdir;
        Tmaxi = intersezione(2)+90;

        if Tmaxi>360
            Tmaxi = Tmaxi - 360;
        end

        Tmini = intersezione(2)-90;

        if Tmini<0
            Tmini = Tmini + 360;
        end

        Condizione_orientazione_reclass = Condizione_orientazione;

        if Tmaxi>0 && Tmaxi<40
            Condizione_orientazione_reclass((Tmaxi<Condizione_orientazione)<Tmini)= NaN;
        else
            Condizione_orientazione_reclass(Condizione_orientazione>Tmaxi)= NaN;
            Condizione_orientazione_reclass(Condizione_orientazione<Tmini)= NaN;
        end

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza>intersezione(1))= NaN;
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        Verifica3D = Condizione_orientazione_reclass & Condizione_pendenza;

        Verifica3D = double(Verifica3D);
        Verifica3Dnonriclass = Verifica3D;


        % write geotiff in output
        filenameout = [fileout3D num2str(Comb_set(j,1)) '-' num2str(Comb_set(j,2))];
        geotiffwrite([DIR_output filenameout],Verifica3D,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VER3D(:,:,j) = Verifica3D;
        VER3Dsum = VER3Dsum + Verifica3Dnonriclass;
    end

    filenameout = [fileout3D 'sum'];
    geotiffwrite([DIR_output filenameout],VER3Dsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% MARKLAND'S TEST FOR TOPPLING

    for i = 1:nset

        % riclassification of aspect
        M_dipdir = Orientazione(i,2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - M_dipdir);
        Tmax = 200;
        Tmin = 160;
                
        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmax)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmin)= NaN;
            
        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Ttop = (90 - Orientazione(i,1))+ Orientazione(i,3);
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<Ttop)= NaN;
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        VerificaTOP = Condizione_orientazione_reclass & Condizione_pendenza;

        VerificaTOP = double(VerificaTOP);
        VerificaTOPnonriclass = VerificaTOP;

        % write geotiff in output
        filenameout = [fileoutTOP num2str(i)];
        geotiffwrite([DIR_output filenameout],VerificaTOP,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VERTOP(:,:,i) = VerificaTOP;
        VERTOPsum = VERTOPsum + VerificaTOPnonriclass;
    end

    filenameout = [fileoutTOP 'sum'];
    geotiffwrite([DIR_output filenameout],VERTOPsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% SUM OF THE RESULTS OF ALL THE POSSIBLE KINEMATISMS
    VERsum = VER2Dsum + VER3Dsum + VERTOPsum;
    geotiffwrite([DIR_output fileoutfinal],VERsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% CREATION OF THE RASTER FILE CONTAINING THE SOURCE POINTS
    VERsource = VERsum;
    VERsource(VERsource>1)= 1;
    geotiffwrite([DIR_output fileoutsource],VERsource,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% end of the process
else
    fprintf('check the input directory: it must contain only aspect and slope files')
end