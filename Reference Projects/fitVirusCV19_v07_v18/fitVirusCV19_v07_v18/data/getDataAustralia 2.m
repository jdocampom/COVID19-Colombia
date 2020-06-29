function [country,C,date0] = getDataAustralia()
%GETDATAAUSTRALIA Coronavirus data for Australia
%  as reported by One World in Data
%     https://ourworldindata.org/coronavirus-source-data
country = 'Australia';
C = [
          6 % 30-Jan-2020
          7 % 31-Jan-2020
          9 % 01-Feb-2020
         11 % 02-Feb-2020
         11 % 03-Feb-2020
         11 % 04-Feb-2020
         12 % 05-Feb-2020
         13 % 06-Feb-2020
         14 % 07-Feb-2020
         14 % 08-Feb-2020
         14 % 09-Feb-2020
         14 % 10-Feb-2020
         14 % 11-Feb-2020
         14 % 12-Feb-2020
         14 % 13-Feb-2020
         15 % 14-Feb-2020
         15 % 15-Feb-2020
         15 % 16-Feb-2020
         15 % 17-Feb-2020
         15 % 18-Feb-2020
         15 % 19-Feb-2020
         15 % 20-Feb-2020
         17 % 21-Feb-2020
         21 % 22-Feb-2020
         22 % 23-Feb-2020
         22 % 24-Feb-2020
         22 % 25-Feb-2020
         22 % 26-Feb-2020
         23 % 27-Feb-2020
         23 % 28-Feb-2020
         25 % 29-Feb-2020
         26 % 01-Mar-2020
         29 % 02-Mar-2020
         33 % 03-Mar-2020
         41 % 04-Mar-2020
         52 % 05-Mar-2020
         59 % 06-Mar-2020
         63 % 07-Mar-2020
         74 % 08-Mar-2020
         80 % 09-Mar-2020
        100 % 10-Mar-2020
        112 % 11-Mar-2020
        126 % 12-Mar-2020
        156 % 13-Mar-2020
        197 % 14-Mar-2020
        249 % 15-Mar-2020
        298 % 16-Mar-2020
        375 % 17-Mar-2020
        454 % 18-Mar-2020
        565 % 19-Mar-2020
        709 % 20-Mar-2020
        874 % 21-Mar-2020
       1098 % 22-Mar-2020
       1709 % 23-Mar-2020
       1823 % 24-Mar-2020
       2423 % 25-Mar-2020
       2799 % 26-Mar-2020
       3166 % 27-Mar-2020
       3378 % 28-Mar-2020
       3809 % 29-Mar-2020
       4093 % 30-Mar-2020
       4557 % 31-Mar-2020
       4707 % 01-Apr-2020
       4976 % 02-Apr-2020
       5224 % 03-Apr-2020
       5548 % 04-Apr-2020
       5687 % 05-Apr-2020
       5744 % 06-Apr-2020
       5844 % 07-Apr-2020
       5956 % 08-Apr-2020
       6052 % 09-Apr-2020
       6152 % 10-Apr-2020
       6238 % 11-Apr-2020
       6289 % 12-Apr-2020
       6322 % 13-Apr-2020
       6366 % 14-Apr-2020
       6416 % 15-Apr-2020
       6458 % 16-Apr-2020
       6497 % 17-Apr-2020
       6533 % 18-Apr-2020
       6586 % 19-Apr-2020
       6612 % 20-Apr-2020
       6625 % 21-Apr-2020
       6647 % 22-Apr-2020
       6654 % 23-Apr-2020
       6667 % 24-Apr-2020
       6687 % 25-Apr-2020
       6703 % 26-Apr-2020
       6713 % 27-Apr-2020
       6725 % 28-Apr-2020
       6738 % 29-Apr-2020
       6746 % 30-Apr-2020
       6762 % 01-May-2020
       6767 % 02-May-2020
       6783 % 03-May-2020
       6801 % 04-May-2020
       6825 % 05-May-2020
       6849 % 06-May-2020
       6875 % 07-May-2020
       6896 % 08-May-2020
       6914 % 09-May-2020
       6929 % 10-May-2020
       6941 % 11-May-2020
%<-------------- add new data here
]';
date0=datenum('30-Jan-2020');
end
