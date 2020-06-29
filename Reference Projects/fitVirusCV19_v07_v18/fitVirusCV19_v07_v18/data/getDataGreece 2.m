function [country,C,date0] = getDataGreece()
%GETDATAGREECE Coronavirus data for Greece
%  as reported by One World in Data
%     https://ourworldindata.org/coronavirus-source-data
country = 'Greece';
C = [
          7 % 01-Mar-2020
          7 % 02-Mar-2020
        NaN % 03-Mar-2020
        NaN % 04-Mar-2020
         10 % 05-Mar-2020
         32 % 06-Mar-2020
         45 % 07-Mar-2020
         66 % 08-Mar-2020
         73 % 09-Mar-2020
         84 % 10-Mar-2020
         90 % 11-Mar-2020
         99 % 12-Mar-2020
        133 % 13-Mar-2020
        190 % 14-Mar-2020
        228 % 15-Mar-2020
        331 % 16-Mar-2020
        352 % 17-Mar-2020
        387 % 18-Mar-2020
        418 % 19-Mar-2020
        464 % 20-Mar-2020
        495 % 21-Mar-2020
        530 % 22-Mar-2020
        624 % 23-Mar-2020
        695 % 24-Mar-2020
        743 % 25-Mar-2020
        821 % 26-Mar-2020
        892 % 27-Mar-2020
        966 % 28-Mar-2020
       1061 % 29-Mar-2020
       1156 % 30-Mar-2020
       1212 % 31-Mar-2020
       1314 % 01-Apr-2020
       1375 % 02-Apr-2020
       1514 % 03-Apr-2020
       1613 % 04-Apr-2020
       1673 % 05-Apr-2020
       1735 % 06-Apr-2020
       1755 % 07-Apr-2020
       1832 % 08-Apr-2020
       1884 % 09-Apr-2020
       1955 % 10-Apr-2020
       2011 % 11-Apr-2020
       2081 % 12-Apr-2020
       2114 % 13-Apr-2020
       2145 % 14-Apr-2020
       2170 % 15-Apr-2020
       2192 % 16-Apr-2020
       2207 % 17-Apr-2020
       2207 % 18-Apr-2020
       2207 % 19-Apr-2020
       2235 % 20-Apr-2020
       2245 % 21-Apr-2020
       2401 % 22-Apr-2020
       2408 % 23-Apr-2020
       2463 % 24-Apr-2020
       2490 % 25-Apr-2020
       2506 % 26-Apr-2020
       2506 % 27-Apr-2020
       2534 % 28-Apr-2020
       2534 % 29-Apr-2020
       2576 % 30-Apr-2020
       2591 % 01-May-2020
       2591 % 02-May-2020
       2620 % 03-May-2020
       2626 % 04-May-2020
       2632 % 05-May-2020
       2642 % 06-May-2020
       2663 % 07-May-2020
       2678 % 08-May-2020
       2691 % 09-May-2020
       2710 % 10-May-2020
       2716 % 11-May-2020
%<-------------- add new data here
]';
date0=datenum('01-Mar-2020');
end
