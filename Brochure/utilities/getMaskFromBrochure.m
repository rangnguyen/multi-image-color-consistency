function mask = getMaskFromBrochure(name)
switch(name)
    case 'Bohol'
        mask = [763        1145          25        1048;
                96         547        1129        2147;
                559         913        1873        2150;
                1015        1693        1132        2151;
                1267        1676        2227        3251;
                1437        2029          28         532;
                99         605        2226        2840;
                708        1164        2229        2914];
    case 'Brazil'
         mask = [1771      2403       268        884;
                1771       2403       986        1608;
                2642       3301       270        886;
                2642       3301       986        1662;
                890        1305       2557       3182;
                890        1305       1830       2374;
                1358       1896       2028       2842;
                2509       3300       1847       3183;
                454        1509       3370       4722;
                2448       3301       3382       4720;];
    case 'Thailand'
         mask = [1         424	      1          493;
                684        991        1          493;
                33         290        521        961;
                297        551        521        707;
                297        551        714        961;
                1          290        989        1233;
                1          290        1239       1485;
                297        515        989        1485;
                749        991        989        1485;
                36         173        1515       1730;
                178        382        1515       1730;
                387        632        1515       1730;
                36         331        1733       1949;
                336        482        1733       1949;
                486        633        1733       1949;];
     case 'Hotel'
         mask = [71        564	      80         786;
                593        1066       526        1224;
                1086       1542       81         783;
                219        600        1238       1811;
                219        600        1829       2407;
                2087       2499       1255       1529;
                2087       2499       1542       1816;
                2087       2499       1828       2101;
                2087       2499       2118       2401;];
    case 'Euro'
         mask = [520       826	      99         581;
                520        826	      603        1083;
                520        826	      1101       1576;
                1064       1371	      99         581;
                1064       1371	      603        1083;
                1064       1371	      1101       1576;
                1748       2053	      99         581;
                1748       2053	      603        1083;
                1748       2053	      1101       1576;];
    case 'Olsen'
         mask = [658       937	      1          280;
                955        1534	      1          280;
                1550       2129	      1          580;
                658        1237	      298        877;
                1253       1532	      596        875;
                1552       1831	      596        1175;
                1850       2129	      596        875;
                1850       2129	      895        1174;
                955        1234	      895        1174;
                1253       1532       895        1474;
                658        937        1192       1471;
                658        937        1490       1769;
                955        1234       1190       1769;
                1253       1532       1490       1769;];
    case 'Birdie_oak'
        mask = [990     1553    1       732;
                990     1553    744     1678;
                93      669     1688    2478;
                677     1553    1686    2478;];
    case 'Appalachian'
        mask = [30      415     1       528;
                807     1209    1       528;
                448     883     541     1042;
                255     978     1072    1564;];
    case 'Surgical'
        mask = [406     611     1       393;
                406     611     401     718;
                406     611     726     1111;
                406     611     1116    1649;];
    case 'Unique_seafood'
        mask = [1123    1517    1279    1701;
                1123    1517    1704    2125;
                1630    2280    1       833;
                1630    2280    836     1647;
                1630    2280    1650    2479;
                1630    2280    2482    3313;
                1630    2280    3316    4129;
                1630    2280    4132    4957;];
    case 'CSH'
        mask = [423     663     1       683;
                686     927     1       683;
                951     1192    1       683;
                135     687     685     1529;
                422     664     1531    1653;
                687     928     1531    1653;
                951     1192    1531    1653;];
    case 'Rhode_island'
        mask = [200     804     1       353;
                942     1542    1       353;
                1438    1701    1068    1313;
                530     819     1625    1874;
                262     1795    1949    2922;];
    case 'FunderMax'
        mask = [120     342     328     604;
                353     543     328     604;
                554     737     328     604;
                120     343     722     957;
                355     605     722     957;
                618     825     722     957;
                120     343     1026    1294;
                355     605     971     1294;
                618     825     1026    1294;
                120     994     1575    2064;];
    case 'Non_surgical'
        mask = [2013    2545    1       829;
                1       638     1127    2197;
                1036    1652    2275    3295;];
    case 'Kashmir'
        mask = [156     444     43      557;
                987     1244    43      285;
                511     836     643     1156;
                102     489     1242    1761;
                661     980     1241    1757;];
    case 'Medical_Assistant'
        mask = [69      468     1       539;
                69      468     541     1058;
                69      468     1060    1597;
                748     1149    970     1573;];
end
        