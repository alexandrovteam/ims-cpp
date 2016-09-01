#include "ms/periodic_table.hpp"

namespace ms {

// values are taken from http://www.ciaaw.org/pubs/TICE-2009.pdf
const std::map<std::string, Element> periodic_table{
    {"H", {"H", 1, {{1.007825032, 2.014101778}, {0.999885, 0.000115}}}},
    {"He", {"He", 2, {{3.016029319, 4.002603254}, {1.34e-06, 0.99999866}}}},
    {"Li", {"Li", 3, {{6.0151223, 7.0160041}, {0.0759, 0.9241}}}},
    {"Be", {"Be", 4, {{9.0121822}, {1.0}}}},
    {"B", {"B", 5, {{10.0129371, 11.0093055}, {0.199, 0.801}}}},
    {"C", {"C", 6, {{12.0, 13.00335484}, {0.9893, 0.0107}}}},
    {"N", {"N", 7, {{14.00307401, 15.00010897}, {0.99636, 0.00364}}}},
    {"O", {"O", 8, {{15.99491462, 16.9991315, 17.9991604}, {0.99757, 0.00038, 0.00205}}}},
    {"F", {"F", 9, {{18.9984032}, {1.0}}}},
    {"Ne",
     {"Ne", 10, {{19.99244018, 20.99384668, 21.99138511}, {0.9048, 0.0027, 0.0925}}}},
    {"Na", {"Na", 11, {{22.98976966}, {1.0}}}},
    {"Mg", {"Mg", 12, {{23.98504187, 24.985837, 25.982593}, {0.7899, 0.1, 0.1101}}}},
    {"Al", {"Al", 13, {{26.98153863}, {1.0}}}},
    {"Si",
     {"Si", 14, {{27.97692649, 28.97649468, 29.97377018}, {0.92223, 0.04685, 0.03092}}}},
    {"P", {"P", 15, {{30.97376149}, {1.0}}}},
    {"S", {"S", 16, {{31.97207073, 32.97145854, 33.96786687, 35.96708088},
                     {0.9499, 0.0075, 0.0425, 0.0001}}}},
    {"Cl", {"Cl", 17, {{34.96885271, 36.9659026}, {0.7576, 0.2424}}}},
    {"Ar", {"Ar", 18,
            {{35.96754511, 37.9627324, 39.96238312}, {0.003336, 0.000629, 0.996035}}}},
    {"K",
     {"K", 19, {{38.9637069, 39.96399867, 40.96182597}, {0.932581, 0.000117, 0.067302}}}},
    {"Ca",
     {"Ca", 20, {{39.9625912, 41.9586183, 42.9587668, 43.9554811, 45.9536927, 47.952533},
                 {0.96941, 0.00647, 0.00135, 0.02086, 4e-05, 0.00187}}}},
    {"Sc", {"Sc", 21, {{44.9559119}, {1.0}}}},
    {"Ti", {"Ti", 22, {{45.9526316, 46.9517631, 47.9479463, 48.94787, 49.9447912},
                       {0.0825, 0.0744, 0.7372, 0.0541, 0.0518}}}},
    {"V", {"V", 23, {{49.9471585, 50.9439595}, {0.0025, 0.9975}}}},
    {"Cr", {"Cr", 24, {{49.9460442, 51.9405075, 52.9406494, 53.9388804},
                       {0.04345, 0.83789, 0.09501, 0.02365}}}},
    {"Mn", {"Mn", 25, {{54.9380451}, {1.0}}}},
    {"Fe", {"Fe", 26, {{53.9396147, 55.9349418, 56.9353983, 57.9332801},
                       {0.05845, 0.91754, 0.02119, 0.00282}}}},
    {"Co", {"Co", 27, {{58.933195}, {1.0}}}},
    {"Ni", {"Ni", 28, {{57.9353429, 59.9307864, 60.931056, 61.9283451, 63.927966},
                       {0.68077, 0.26223, 0.011399, 0.036346, 0.009255}}}},
    {"Cu", {"Cu", 29, {{62.9295975, 64.9277895}, {0.6915, 0.3085}}}},
    {"Zn", {"Zn", 30, {{63.9291422, 65.9260334, 66.9271273, 67.9248442, 69.9253193},
                       {0.4917, 0.2773, 0.0404, 0.1845, 0.0061}}}},
    {"Ga", {"Ga", 31, {{68.9255736, 70.9247013}, {0.60108, 0.39892}}}},
    {"Ge", {"Ge", 32, {{69.9242474, 71.9220758, 72.9234589, 73.9211778, 75.9214026},
                       {0.2057, 0.2745, 0.0775, 0.3650, 0.0773}}}},
    {"As", {"As", 33, {{74.9215965}, {1.0}}}},
    {"Se",
     {"Se", 34, {{73.9224764, 75.9192136, 76.919914, 77.9173091, 79.9165213, 81.9166994},
                 {0.0089, 0.0937, 0.0763, 0.2377, 0.4961, 0.0873}}}},
    {"Br", {"Br", 35, {{78.9183379, 80.916291}, {0.5069, 0.4931}}}},
    {"Kr",
     {"Kr", 36, {{77.9203648, 79.916379, 81.9134836, 82.914136, 83.911507, 85.91061073},
                 {0.00355, 0.02286, 0.11593, 0.115, 0.56987, 0.17279}}}},
    {"Rb", {"Rb", 37, {{84.91178974, 86.90918053}, {0.7217, 0.2783}}}},
    {"Sr", {"Sr", 38, {{83.913425, 85.9092602, 86.9088771, 87.9056121},
                       {0.0056, 0.0986, 0.07, 0.8258}}}},
    {"Y", {"Y", 39, {{88.9058483}, {1.0}}}},
    {"Zr", {"Zr", 40, {{89.9047044, 90.9056458, 91.9050408, 93.9063152, 95.9082734},
                       {0.5145, 0.1122, 0.1715, 0.1738, 0.028}}}},
    {"Nb", {"Nb", 41, {{92.9063781}, {1.0}}}},
    {"Mo", {"Mo", 42, {{91.906811, 93.9050883, 94.9058421, 95.9046795, 96.9060215,
                        97.9054082, 99.90747},
                       {0.1453, 0.0915, 0.1584, 0.1667, 0.0960, 0.2439, 0.0982}}}},
    {"Tc", {"Tc", 43, {{96.9064}, {1.0}}}},
    {"Ru", {"Ru", 44, {{95.907598, 97.905287, 98.9059393, 99.9042195, 100.9055821,
                        101.9043493, 103.905433},
                       {0.0554, 0.0187, 0.1276, 0.126, 0.1706, 0.3155, 0.1862}}}},
    {"Rh", {"Rh", 45, {{102.905504}, {1.0}}}},
    {"Pd",
     {"Pd", 46, {{101.905609, 103.904036, 104.905085, 105.903486, 107.903892, 109.905153},
                 {0.0102, 0.1114, 0.2233, 0.2733, 0.2646, 0.1172}}}},
    {"Ag", {"Ag", 47, {{106.905097, 108.904752}, {0.51839, 0.48161}}}},
    {"Cd", {"Cd", 48, {{105.906459, 107.904184, 109.9030021, 110.9041781, 111.9027578,
                        112.9044017, 113.9033585, 115.904756},
                       {0.0125, 0.0089, 0.1249, 0.128, 0.2413, 0.1222, 0.2873, 0.0749}}}},
    {"In", {"In", 49, {{112.904058, 114.903878}, {0.0429, 0.9571}}}},
    {"Sn", {"Sn", 50, {{111.904818, 113.902779, 114.903342, 115.901741, 116.902952,
                        117.901603, 118.903308, 119.9021947, 121.903439, 123.9052739},
                       {0.0097, 0.0066, 0.0034, 0.1454, 0.0768, 0.2422, 0.0859, 0.3258,
                        0.0463, 0.0579}}}},
    {"Sb", {"Sb", 51, {{120.9038157, 122.904214}, {0.5721, 0.4279}}}},
    {"Te",
     {"Te", 52, {{119.90402, 121.9030439, 122.90427, 123.9028179, 124.9044307,
                  125.9033117, 127.9044631, 129.9062244},
                 {0.0009, 0.0255, 0.0089, 0.0474, 0.0707, 0.1884, 0.3174, 0.3408}}}},
    {"I", {"I", 53, {{126.904473}, {1.0}}}},
    {"Xe", {"Xe", 54, {{123.905893, 125.904274, 127.9035313, 128.9047794, 129.903508,
                        130.9050824, 131.9041535, 133.9053945, 135.907219},
                       {0.000952, 0.00089, 0.019102, 0.264006, 0.04071, 0.212324,
                        0.269086, 0.104357, 0.088573}}}},
    {"Cs", {"Cs", 55, {{132.9054519}, {1.0}}}},
    {"Ba", {"Ba", 56, {{129.9063208, 131.9050613, 133.9045084, 134.9056886, 135.9045759,
                        136.9058274, 137.9052472},
                       {0.00106, 0.00101, 0.02417, 0.06592, 0.07854, 0.11232, 0.71698}}}},
    {"La", {"La", 57, {{137.907112, 138.9063533}, {0.0008881, 0.9991119}}}},
    {"Ce", {"Ce", 58, {{135.907172, 137.905991, 139.9054387, 141.909244},
                       {0.00185, 0.00251, 0.8845, 0.11114}}}},
    {"Pr", {"Pr", 59, {{140.9076528}, {1.0}}}},
    {"Nd", {"Nd", 60, {{141.9077233, 142.9098143, 143.9100873, 144.9125736, 145.9131169,
                        147.916893, 149.920891},
                       {0.27152, 0.12174, 0.23798, 0.08293, 0.17189, 0.05756, 0.05638}}}},
    {"Pm", {"Pm", 61, {{144.9127}, {1.0}}}},
    {"Sm", {"Sm", 62, {{143.911999, 146.9148979, 147.9148227, 148.9171847, 149.9172755,
                        151.9197324, 153.9222093},
                       {0.0307, 0.1499, 0.1124, 0.1382, 0.0738, 0.2675, 0.2275}}}},
    {"Eu", {"Eu", 63, {{150.9198502, 152.9212303}, {0.4781, 0.5219}}}},
    {"Gd", {"Gd", 64, {{151.919791, 153.9208656, 154.922622, 155.9221227, 156.9239601,
                        157.9241039, 159.9270541},
                       {0.002, 0.0218, 0.148, 0.2047, 0.1565, 0.2484, 0.2186}}}},
    {"Tb", {"Tb", 65, {{158.9253468}, {1.0}}}},
    {"Dy", {"Dy", 66, {{155.924283, 157.924409, 159.9251975, 160.9269334, 161.9267984,
                        162.9287312, 163.9291748},
                       {0.00056, 0.00095, 0.02329, 0.18889, 0.25475, 0.24896, 0.2826}}}},
    {"Ho", {"Ho", 67, {{164.9303221}, {1.0}}}},
    {"Er", {"Er", 68,
            {{161.928778, 163.9292, 165.9302931, 166.9320482, 167.9323702, 169.9354643},
             {0.00139, 0.01601, 0.33503, 0.22869, 0.26978, 0.1491}}}},
    {"Tm", {"Tm", 69, {{168.9342133}, {1.0}}}},
    {"Yb", {"Yb", 70, {{167.933897, 169.9347618, 170.9363258, 171.9363815, 172.9382108,
                        173.9388621, 175.9425717},
                       {0.00123, 0.02982, 0.1409, 0.2168, 0.16103, 0.32026, 0.12996}}}},
    {"Lu", {"Lu", 71, {{174.9407718, 175.9426863}, {0.97401, 0.02599}}}},
    {"Hf", {"Hf", 72,
            {{173.940046, 175.9414086, 176.9432207, 177.9436988, 178.9458161, 179.94655},
             {0.0016, 0.0526, 0.186, 0.2728, 0.1362, 0.3508}}}},
    {"Ta", {"Ta", 73, {{179.9474648, 180.9479958}, {0.0001201, 0.9998799}}}},
    {"W", {"W", 74, {{179.946704, 181.9482042, 182.950223, 183.9509312, 185.9543641},
                     {0.0012, 0.265, 0.1431, 0.3064, 0.2843}}}},
    {"Re", {"Re", 75, {{184.952955, 186.9557531}, {0.374, 0.626}}}},
    {"Os", {"Os", 76, {{183.9524891, 185.9538382, 186.9557505, 187.9558382, 188.9581475,
                        189.958447, 191.9614807},
                       {0.0002, 0.0159, 0.0196, 0.1324, 0.1615, 0.2626, 0.4078}}}},
    {"Ir", {"Ir", 77, {{190.960594, 192.9629264}, {0.373, 0.627}}}},
    {"Pt", {"Pt", 78,
            {{189.959932, 191.961038, 193.9626803, 194.9647911, 195.9649515, 197.967893},
             {0.00012, 0.00782, 0.3286, 0.3378, 0.2521, 0.07356}}}},
    {"Au", {"Au", 79, {{196.9665687}, {1.0}}}},
    {"Hg", {"Hg", 80, {{195.965833, 197.966769, 198.9682799, 199.968326, 200.9703023,
                        201.970643, 203.9734939},
                       {0.0015, 0.0997, 0.1687, 0.231, 0.1318, 0.2986, 0.0687}}}},
    {"Tl", {"Tl", 81, {{202.9723442, 204.9744275}, {0.2952, 0.7048}}}},
    {"Pb", {"Pb", 82, {{203.9730436, 205.9744653, 206.9758969, 207.9766521},
                       {0.014, 0.241, 0.221, 0.524}}}},
    {"Bi", {"Bi", 83, {{208.9803987}, {1.0}}}}, {"Po", {"Po", 84, {{209.0}, {1.0}}}},
    {"At", {"At", 85, {{210.0}, {1.0}}}}, {"Rn", {"Rn", 86, {{220.0}, {1.0}}}},
    {"Fr", {"Fr", 87, {{223.0}, {1.0}}}}, {"Ra", {"Ra", 88, {{226.0}, {1.0}}}},
    {"Ac", {"Ac", 89, {{227.0}, {1.0}}}}, {"Th", {"Th", 90, {{232.0380553}, {1.0}}}},
    {"Pa", {"Pa", 91, {{231.035884}, {1.0}}}},
    {"U",
     {"U", 92, {{234.0409521, 235.0439299, 238.0507882}, {5.4e-05, 0.007204, 0.992742}}}},
    {"Np", {"Np", 93, {{237.0}, {1.0}}}}, {"Pu", {"Pu", 94, {{244.0}, {1.0}}}},
    {"Am", {"Am", 95, {{243.0}, {1.0}}}}, {"Cm", {"Cm", 96, {{247.0}, {1.0}}}},
    {"Bk", {"Bk", 97, {{247.0}, {1.0}}}}, {"Cf", {"Cf", 98, {{251.0}, {1.0}}}},
    {"Es", {"Es", 99, {{252.0}, {1.0}}}}, {"Fm", {"Fm", 100, {{257.0}, {1.0}}}},
    {"Md", {"Md", 101, {{258.0}, {1.0}}}}, {"No", {"No", 102, {{259.0}, {1.0}}}},
    {"Lr", {"Lr", 103, {{262.0}, {1.0}}}}, {"Rf", {"Rf", 104, {{261.0}, {1.0}}}},
    {"Db", {"Db", 105, {{262.0}, {1.0}}}}, {"Sg", {"Sg", 106, {{266.0}, {1.0}}}},
};
}
