function [h]=plotCIELabMap(L_star,a_star_axis,b_star_axis,varargin)
%h=plotCIELabMap(L_star,a_star,b_star,options)
%
%Produces a uniform lightness colormap given ranges of a and be coordiantes
%using sRGB values
%
%Input values
%   lightness value     L_star=80;
%   a-axis               a_star=[-40:40];
%   b-axis               b_star=[-40:40];
%
%Example:
%
%plotCIELabMap(60,[-40:2:40],[-40:2:40])
%
%Options 
%Name-Value Arguments
%         ColorSpace - Color space of the output RGB values
%             'srgb' (default) | 'adobe-rgb-1998' | 'prophoto-rgb' |
%             'linear-rgb'
%         WhitePoint - Reference white point
%             'a' | 'c' | 'e' | 'd50' | 'd55' | 'd65' (default) | 'icc' |
%             1-by-3 vector
%
%The function saves the last generated map to the userpath()
%
%   To redo the calculation run with 'reset' as an argument: Example:
%
% plotCIELabMap(80,[-40:40],[-40:40], 'ColorSpace', 'srgb','WhitePoint','d65','reset')
%
%
% Anders Thorseth, DTU Electro, 2025


% Define optional parameters for lab2rgb conversion
p = inputParser;
addParameter(p, 'ColorSpace', 'srgb');   % Default value = 'srgb'
addParameter(p, 'WhitePoint','d65');    % Default value = 'd65'


filePathCIELabColorMap=fullfile(userpath,'CIELabColorMap'); %path were the data is stored for later use

persistent Lab_array
recalculate=false; % initially value

vararginNew=varargin;

%use 'reset' as an agument to clean up
if sum(strcmp(varargin,'reset'))>0
    Lab_array=[];
    delete(filePathCIELabColorMap)
    disp('deleting old data')
    idx = strcmp(vararginNew, 'reset');
    vararginNew(idx) = [];
end

%triggers recalcualtion or load if 'Lab_array' is empty
if isempty(Lab_array)
    recalculate=true;
    if exist(filePathCIELabColorMap,"file")
        disp('loading from file')
        load(filePathCIELabColorMap,"Lab_array","a_star_axis","b_star_axis"); % load data for Lab_array
        recalculate=false;
    end
end

if recalculate

    disp('recalculating')

    %set default values 
    if nargin<1
        L_star=100;
    end

    if nargin<3
        a_star_axis=-40:40;
        b_star_axis=-40:40;
    end

    [a_grid,b_grid]=meshgrid(a_star_axis,b_star_axis); %set up a,b grid points

    % WhitePoint — Reference white point
    % "d65" (default) | "a" | "c" | "e" | "d50" | "d55" | "icc" | 1-by-3 vector
    Lab_array=zeros(length(a_grid),length(a_grid),3); %initiate matrix
    for n=1:length(a_star_axis)
        for m=1:length(b_star_axis)
            CIE_Lab=[L_star,a_grid(n,m),b_grid(n,m)];
            Lab_array(n,m,:) =lab2rgb(CIE_Lab,vararginNew{:});
        end
    end
    save(fullfile(userpath,'CIELabColorMap'),"Lab_array","a_star_axis","b_star_axis")
end


h=image(a_star_axis,b_star_axis,Lab_array);
xlabel("CIE a*")
ylabel("CIE b*")
set(gca,"YDir",'normal')
