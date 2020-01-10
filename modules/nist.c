#include <string.h>
#include <stdlib.h>

#include "globals.h"
#include "nist.h"

/******************************************************************************

 Function nist_atomic_mass(): return the atomic mass (in amu) for a given
 isotope. See the isotope enumeration in the header file for details.

 Date: 14/10/2017
 Font: NIST Physical Measurement Laboratory
 Link: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all

******************************************************************************/

double nist_atomic_mass(const isotope a)
{
	switch (a)
	{
		case    0: return 1.007825032239;   // 1H
		case    1: return 2.0141017781212;  // 2H
		case    2: return 3.016049277924;   // 3H
		case    3: return 4.0264311;        // 4H
		case    4: return 5.03531196;       // 5H
		case    5: return 6.0449627;        // 6H
		case    6: return 7.052711;         // 7H
		case    7: return 3.016029320125;   // 3He
		case    8: return 4.002603254136;   // 4He
		case    9: return 5.01205721;       // 5He
		case   10: return 6.01888589157;    // 6He
		case   11: return 7.027990781;      // 7He
		case   12: return 8.03393439095;    // 8He
		case   13: return 9.04394650;       // 9He
		case   14: return 10.0527911;       // 10He
		case   15: return 3.030821;         // 3Li
		case   16: return 4.0271923;        // 4Li
		case   17: return 5.01253854;       // 5Li
		case   18: return 6.015122887416;   // 6Li
		case   19: return 7.016003436645;   // 7Li
		case   20: return 8.02248624650;    // 8Li
		case   21: return 9.0267901920;     // 9Li
		case   22: return 10.03548314;      // 10Li
		case   23: return 11.0437235866;    // 11Li
		case   24: return 12.05251716;      // 12Li
		case   25: return 13.0626338;       // 13Li
		case   26: return 5.039922;         // 5Be
		case   27: return 6.019726458;      // 6Be
		case   28: return 7.01692871776;    // 7Be
		case   29: return 8.00530510237;    // 8Be
		case   30: return 9.01218306582;    // 9Be
		case   31: return 10.01353469586;   // 10Be
		case   32: return 11.0216610826;    // 11Be
		case   33: return 12.026922120;     // 12Be
		case   34: return 13.03613511;      // 13Be
		case   35: return 14.0428914;       // 14Be
		case   36: return 15.0534243;       // 15Be
		case   37: return 16.0616718;       // 16Be
		case   38: return 6.050822;         // 6B
		case   39: return 7.02971227;       // 7B
		case   40: return 8.024607311;      // 8B
		case   41: return 9.0133296597;     // 9B
		case   42: return 10.0129369541;    // 10B
		case   43: return 11.0093053645;    // 11B
		case   44: return 12.014352714;     // 12B
		case   45: return 13.017780212;     // 13B
		case   46: return 14.02540423;      // 14B
		case   47: return 15.03108823;      // 15B
		case   48: return 16.03984226;      // 16B
		case   49: return 17.0469918;       // 17B
		case   50: return 18.0556618;       // 18B
		case   51: return 19.0631043;       // 19B
		case   52: return 20.0720775;       // 20B
		case   53: return 21.0812997;       // 21B
		case   54: return 8.03764320;       // 8C
		case   55: return 9.031037223;      // 9C
		case   56: return 10.0168533142;    // 10C
		case   57: return 11.011433610;     // 11C
		case   58: return 12.000000000;     // 12C
		case   59: return 13.0033548350723; // 13C
		case   60: return 14.003241988440;  // 14C
		case   61: return 15.0105992686;    // 15C
		case   62: return 16.014701338;     // 16C
		case   63: return 17.02257719;      // 17C
		case   64: return 18.02675132;      // 18C
		case   65: return 19.0348011;       // 19C
		case   66: return 20.0403226;       // 20C
		case   67: return 21.0490043;       // 21C
		case   68: return 22.0575326;       // 22C
		case   69: return 23.068911;        // 23C
		case   70: return 10.0416543;       // 10N
		case   71: return 11.02609150;      // 11N
		case   72: return 12.018613211;     // 12N
		case   73: return 13.0057386129;    // 13N
		case   74: return 14.0030740044320; // 14N
		case   75: return 15.0001088988864; // 15N
		case   76: return 16.006101925;     // 16N
		case   77: return 17.00844916;      // 17N
		case   78: return 18.01407820;      // 18N
		case   79: return 19.01702218;      // 19N
		case   80: return 20.02336660;      // 20N
		case   81: return 21.0271110;       // 21N
		case   82: return 22.0343921;       // 22N
		case   83: return 23.0411432;       // 23N
		case   84: return 24.0503943;       // 24N
		case   85: return 25.0601054;       // 25N
		case   86: return 12.03426226;      // 12O
		case   87: return 13.02481510;      // 13O
		case   88: return 14.0085963612;    // 14O
		case   89: return 15.0030656253;    // 15O
		case   90: return 15.9949146195717; // 16O
		case   91: return 16.9991317565069; // 17O
		case   92: return 17.9991596128676; // 18O
		case   93: return 19.003578028;     // 19O
		case   94: return 20.0040753595;    // 20O
		case   95: return 21.00865513;      // 21O
		case   96: return 22.00996661;      // 22O
		case   97: return 23.01569697;      // 23O
		case   98: return 24.0198612;       // 24O
		case   99: return 25.0293612;       // 25O
		case  100: return 26.0372917;       // 26O
		case  101: return 27.0477254;       // 27O
		case  102: return 28.0559175;       // 28O
		case  103: return 14.03431544;      // 14F
		case  104: return 15.01804367;      // 15F
		case  105: return 16.011465789;     // 16F
		case  106: return 17.0020952427;    // 17F
		case  107: return 18.0009373350;    // 18F
		case  108: return 18.9984031627392; // 19F
		case  109: return 19.99998125231;   // 20F
		case  110: return 20.999948919;     // 21F
		case  111: return 22.00299913;      // 22F
		case  112: return 23.00355754;      // 23F
		case  113: return 24.00811578;      // 24F
		case  114: return 25.01219981;      // 25F
		case  115: return 26.02003883;      // 26F
		case  116: return 27.0264420;       // 27F
		case  117: return 28.0353421;       // 28F
		case  118: return 29.0425454;       // 29F
		case  119: return 30.0516564;       // 30F
		case  120: return 31.0597156;       // 31F
		case  121: return 16.02575022;      // 16Ne
		case  122: return 17.0177139638;    // 17Ne
		case  123: return 18.0057087039;    // 18Ne
		case  124: return 19.0018809117;    // 19Ne
		case  125: return 19.992440176217;  // 20Ne
		case  126: return 20.99384668541;   // 21Ne
		case  127: return 21.99138511418;   // 22Ne
		case  128: return 22.9944669111;    // 23Ne
		case  129: return 23.9936106555;    // 24Ne
		case  130: return 24.99778948;      // 25Ne
		case  131: return 26.00051520;      // 26Ne
		case  132: return 27.00755370;      // 27Ne
		case  133: return 28.0121210;       // 28Ne
		case  134: return 29.0197511;       // 29Ne
		case  135: return 30.0247330;       // 30Ne
		case  136: return 31.033117;        // 31Ne
		case  137: return 32.0397254;       // 32Ne
		case  138: return 33.0493864;       // 33Ne
		case  139: return 34.0567355;       // 34Ne
		case  140: return 18.0268812;       // 18Na
		case  141: return 19.01388011;      // 19Na
		case  142: return 20.007354412;     // 20Na
		case  143: return 20.9976546930;    // 21Na
		case  144: return 21.9944374118;    // 22Na
		case  145: return 22.989769282019;  // 23Na
		case  146: return 23.99096295038;   // 24Na
		case  147: return 24.989954013;     // 25Na
		case  148: return 25.992634638;     // 26Na
		case  149: return 26.994076540;     // 27Na
		case  150: return 27.99893911;      // 28Na
		case  151: return 29.002877179;     // 29Na
		case  152: return 30.009097951;     // 30Na
		case  153: return 31.01316325;      // 31Na
		case  154: return 32.0201913;       // 32Na
		case  155: return 33.0257364;       // 33Na
		case  156: return 34.0335954;       // 34Na
		case  157: return 35.0406263;       // 35Na
		case  158: return 36.0492964;       // 36Na
		case  159: return 37.0570565;       // 37Na
		case  160: return 19.03416954;      // 19Mg
		case  161: return 20.01885029;      // 20Mg
		case  162: return 21.01171618;      // 21Mg
		case  163: return 21.9995706534;    // 22Mg
		case  164: return 22.9941242174;    // 23Mg
		case  165: return 23.98504169714;   // 24Mg
		case  166: return 24.98583697650;   // 25Mg
		case  167: return 25.98259296831;   // 26Mg
		case  168: return 26.98434062453;   // 27Mg
		case  169: return 27.983876722;     // 28Mg
		case  170: return 28.98861712;      // 29Mg
		case  171: return 29.990462937;     // 30Mg
		case  172: return 30.996648033;     // 31Mg
		case  173: return 31.999110234;     // 32Mg
		case  174: return 33.005327131;     // 33Mg
		case  175: return 34.00893531;      // 34Mg
		case  176: return 35.0167919;       // 35Mg
		case  177: return 36.0218849;       // 36Mg
		case  178: return 37.0303754;       // 37Mg
		case  179: return 38.0365854;       // 38Mg
		case  180: return 39.0453855;       // 39Mg
		case  181: return 40.0521864;       // 40Mg
		case  182: return 21.0289743;       // 21Al
		case  183: return 22.0195443;       // 22Al
		case  184: return 23.0072443537;    // 23Al
		case  185: return 23.999948912;     // 24Al
		case  186: return 24.9904281051;    // 25Al
		case  187: return 25.98689190469;   // 26Al
		case  188: return 26.9815385311;    // 27Al
		case  189: return 27.9819102113;    // 28Al
		case  190: return 28.980456510;     // 29Al
		case  191: return 29.98296015;      // 30Al
		case  192: return 30.98394522;      // 31Al
		case  193: return 31.98808513;      // 32Al
		case  194: return 32.99090981;      // 33Al
		case  195: return 33.99670574;      // 34Al
		case  196: return 34.99976475;      // 35Al
		case  197: return 36.0063911;       // 36Al
		case  198: return 37.0105313;       // 37Al
		case  199: return 38.0174027;       // 38Al
		case  200: return 39.0225454;       // 39Al
		case  201: return 40.0300354;       // 40Al
		case  202: return 41.0363864;       // 41Al
		case  203: return 42.0438464;       // 42Al
		case  204: return 43.0514775;       // 43Al
		case  205: return 22.0357954;       // 22Si
		case  206: return 23.0254454;       // 23Si
		case  207: return 24.01153521;      // 24Si
		case  208: return 25.00410911;      // 25Si
		case  209: return 25.9923338411;    // 26Si
		case  210: return 26.9867048115;    // 27Si
		case  211: return 27.9769265346544; // 28Si
		case  212: return 28.9764946649052; // 29Si
		case  213: return 29.97377013623;   // 30Si
		case  214: return 30.97536319446;   // 31Si
		case  215: return 31.9741515432;    // 32Si
		case  216: return 32.9779769675;    // 33Si
		case  217: return 33.97857615;      // 34Si
		case  218: return 34.98458341;      // 35Si
		case  219: return 35.98669577;      // 36Si
		case  220: return 36.99292189;      // 37Si
		case  221: return 37.99552375;      // 38Si
		case  222: return 39.00249197;      // 39Si
		case  223: return 40.0058325;       // 40Si
		case  224: return 41.0130140;       // 41Si
		case  225: return 42.0177854;       // 42Si
		case  226: return 43.0248064;       // 43Si
		case  227: return 44.0306164;       // 44Si
		case  228: return 45.0399575;       // 45Si
		case  229: return 24.0357754;       // 24P
		case  230: return 25.0211943;       // 25P
		case  231: return 26.0117821;       // 26P
		case  232: return 26.99922428;      // 27P
		case  233: return 27.992326612;     // 28P
		case  234: return 28.9818007960;    // 29P
		case  235: return 29.9783137534;    // 30P
		case  236: return 30.9737619984270; // 31P
		case  237: return 31.97390764342;   // 32P
		case  238: return 32.971725712;     // 33P
		case  239: return 33.9736458987;    // 34P
		case  240: return 34.973314120;     // 35P
		case  241: return 35.97826014;      // 36P
		case  242: return 36.97960741;      // 37P
		case  243: return 37.98425293;      // 38P
		case  244: return 38.98622798;      // 39P
		case  245: return 39.9913312;       // 40P
		case  246: return 40.99465486;      // 41P
		case  247: return 42.0010823;       // 42P
		case  248: return 43.0050240;       // 43P
		case  249: return 44.0112154;       // 44P
		case  250: return 45.0164564;       // 45P
		case  251: return 46.0244675;       // 46P
		case  252: return 47.0313986;       // 47P
		case  253: return 26.0290764;       // 26S
		case  254: return 27.0182843;       // 27S
		case  255: return 28.0043717;       // 28S
		case  256: return 28.99661154;      // 29S
		case  257: return 29.9849070340;    // 30S
		case  258: return 30.9795570125;    // 31S
		case  259: return 31.972071174414;  // 32S
		case  260: return 32.971458909815;  // 33S
		case  261: return 33.96786700447;   // 34S
		case  262: return 34.96903231043;   // 35S
		case  263: return 35.9670807120;    // 36S
		case  264: return 36.9711255121;    // 37S
		case  265: return 37.971163377;     // 38S
		case  266: return 38.97513454;      // 39S
		case  267: return 39.975482643;     // 40S
		case  268: return 40.979593544;     // 41S
		case  269: return 41.981065130;     // 42S
		case  270: return 42.986907653;     // 43S
		case  271: return 43.990118856;     // 44S
		case  272: return 44.9957274;       // 45S
		case  273: return 46.0000454;       // 46S
		case  274: return 47.0079554;       // 47S
		case  275: return 48.0137064;       // 48S
		case  276: return 49.0227672;       // 49S
		case  277: return 28.0295464;       // 28Cl
		case  278: return 29.0147843;       // 29Cl
		case  279: return 30.0047721;       // 30Cl
		case  280: return 30.99241454;      // 31Cl
		case  281: return 31.9856846460;    // 32Cl
		case  282: return 32.9774519942;    // 33Cl
		case  283: return 33.97376248552;   // 34Cl
		case  284: return 34.96885268237;   // 35Cl
		case  285: return 35.96830680938;   // 36Cl
		case  286: return 36.96590260255;   // 37Cl
		case  287: return 37.9680104411;    // 38Cl
		case  288: return 38.968008219;     // 39Cl
		case  289: return 39.97041534;      // 40Cl
		case  290: return 40.97068574;      // 41Cl
		case  291: return 41.9732515;       // 42Cl
		case  292: return 42.9738910;       // 43Cl
		case  293: return 43.9778720;       // 44Cl
		case  294: return 44.9802911;       // 45Cl
		case  295: return 45.9851717;       // 46Cl
		case  296: return 46.9891643;       // 47Cl
		case  297: return 47.9956454;       // 48Cl
		case  298: return 49.0012364;       // 49Cl
		case  299: return 50.0090564;       // 50Cl
		case  300: return 51.0155475;       // 51Cl
		case  301: return 30.0230754;       // 30Ar
		case  302: return 31.0121222;       // 31Ar
		case  303: return 31.997637819;     // 32Ar
		case  304: return 32.9899255543;    // 33Ar
		case  305: return 33.98027009083;   // 34Ar
		case  306: return 34.9752575980;    // 35Ar
		case  307: return 35.96754510528;   // 36Ar
		case  308: return 36.9667763322;    // 37Ar
		case  309: return 37.9627321121;    // 38Ar
		case  310: return 38.964313054;     // 39Ar
		case  311: return 39.962383123724;  // 40Ar
		case  312: return 40.9645005737;    // 41Ar
		case  313: return 41.963045762;     // 42Ar
		case  314: return 42.965636157;     // 43Ar
		case  315: return 43.964923817;     // 44Ar
		case  316: return 44.9680397355;    // 45Ar
		case  317: return 45.96808344;      // 46Ar
		case  318: return 46.97293596;      // 47Ar
		case  319: return 47.9759132;       // 48Ar
		case  320: return 48.9819043;       // 49Ar
		case  321: return 49.9861354;       // 50Ar
		case  322: return 50.9937064;       // 51Ar
		case  323: return 51.9989664;       // 52Ar
		case  324: return 53.0072975;       // 53Ar
		case  325: return 32.0226554;       // 32K
		case  326: return 33.0075621;       // 33K
		case  327: return 33.9986932;       // 34K
		case  328: return 34.9880054155;    // 35K
		case  329: return 35.9813020137;    // 36K
		case  330: return 36.9733758910;    // 37K
		case  331: return 37.9690811221;    // 38K
		case  332: return 38.963706486449;  // 39K
		case  333: return 39.96399816660;   // 40K
		case  334: return 40.961825257941;  // 41K
		case  335: return 41.9624023111;    // 42K
		case  336: return 42.9607347044;    // 43K
		case  337: return 43.9615869945;    // 44K
		case  338: return 44.9606914956;    // 45K
		case  339: return 45.9619815978;    // 46K
		case  340: return 46.961661615;     // 47K
		case  341: return 47.9653411983;    // 48K
		case  342: return 48.9682107586;    // 49K
		case  343: return 49.972380083;     // 50K
		case  344: return 50.97582814;      // 51K
		case  345: return 51.9822443;       // 52K
		case  346: return 52.9874654;       // 53K
		case  347: return 53.9946364;       // 54K
		case  348: return 55.0007675;       // 55K
		case  349: return 56.0085186;       // 56K
		case  350: return 34.0148732;       // 34Ca
		case  351: return 35.0051421;       // 35Ca
		case  352: return 35.99307443;      // 36Ca
		case  353: return 36.9858978568;    // 37Ca
		case  354: return 37.9763192221;    // 38Ca
		case  355: return 38.9707108164;    // 39Ca
		case  356: return 39.96259086322;   // 40Ca
		case  357: return 40.9622779215;    // 41Ca
		case  358: return 41.9586178316;    // 42Ca
		case  359: return 42.9587664424;    // 43Ca
		case  360: return 43.9554815635;    // 44Ca
		case  361: return 44.9561863539;    // 45Ca
		case  362: return 45.953689024;     // 46Ca
		case  363: return 46.954542424;     // 47Ca
		case  364: return 47.9525227613;    // 48Ca
		case  365: return 48.9556627423;    // 49Ca
		case  366: return 49.957499217;     // 50Ca
		case  367: return 50.96098924;      // 51Ca
		case  368: return 51.96321764;      // 52Ca
		case  369: return 52.9694543;       // 53Ca
		case  370: return 53.9734054;       // 54Ca
		case  371: return 54.9803054;       // 55Ca
		case  372: return 55.9850864;       // 56Ca
		case  373: return 56.9926264;       // 57Ca
		case  374: return 57.9979475;       // 58Ca
		case  375: return 36.0164832;       // 36Sc
		case  376: return 37.0037432;       // 37Sc
		case  377: return 37.9951221;       // 38Sc
		case  378: return 38.98478526;      // 39Sc
		case  379: return 39.977967330;     // 40Sc
		case  380: return 40.96925110588;   // 41Sc
		case  381: return 41.9655165318;    // 42Sc
		case  382: return 42.961150520;     // 43Sc
		case  383: return 43.959402919;     // 44Sc
		case  384: return 44.9559082877;    // 45Sc
		case  385: return 45.9551682678;    // 46Sc
		case  386: return 46.952403721;     // 47Sc
		case  387: return 47.952223653;     // 48Sc
		case  388: return 48.950014629;     // 49Sc
		case  389: return 49.95217616;      // 50Sc
		case  390: return 50.95359221;      // 51Sc
		case  391: return 51.9568815;       // 52Sc
		case  392: return 52.9590929;       // 53Sc
		case  393: return 53.9639339;       // 54Sc
		case  394: return 54.9678250;       // 55Sc
		case  395: return 55.9734543;       // 56Sc
		case  396: return 56.9777754;       // 57Sc
		case  397: return 57.9840364;       // 58Sc
		case  398: return 58.9889464;       // 59Sc
		case  399: return 59.9956575;       // 60Sc
		case  400: return 61.0010086;       // 61Sc
		case  401: return 38.0114532;       // 38Ti
		case  402: return 39.0023622;       // 39Ti
		case  403: return 39.9905017;       // 40Ti
		case  404: return 40.98314830;      // 41Ti
		case  405: return 41.9730490330;    // 42Ti
		case  406: return 42.968522578;     // 43Ti
		case  407: return 43.9596899575;    // 44Ti
		case  408: return 44.9581219895;    // 45Ti
		case  409: return 45.9526277235;    // 46Ti
		case  410: return 46.9517587938;    // 47Ti
		case  411: return 47.9479419838;    // 48Ti
		case  412: return 48.9478656839;    // 49Ti
		case  413: return 49.9447868939;    // 50Ti
		case  414: return 50.9466106565;    // 51Ti
		case  415: return 51.946893076;     // 52Ti
		case  416: return 52.9497311;       // 53Ti
		case  417: return 53.9510513;       // 54Ti
		case  418: return 54.9552717;       // 55Ti
		case  419: return 55.9579115;       // 56Ti
		case  420: return 56.9636427;       // 57Ti
		case  421: return 57.9666043;       // 58Ti
		case  422: return 58.9724743;       // 59Ti
		case  423: return 59.9760354;       // 60Ti
		case  424: return 60.9824564;       // 61Ti
		case  425: return 61.9865175;       // 62Ti
		case  426: return 62.9937575;       // 63Ti
		case  427: return 40.0127643;       // 40V
		case  428: return 41.0002132;       // 41V
		case  429: return 41.9918232;       // 42V
		case  430: return 42.98076646;      // 43V
		case  431: return 43.9741120;       // 44V
		case  432: return 44.965774886;     // 45V
		case  433: return 45.9601987836;    // 46V
		case  434: return 46.9549049136;    // 47V
		case  435: return 47.952252211;     // 48V
		case  436: return 48.9485118096;    // 49V
		case  437: return 49.9471560195;    // 50V
		case  438: return 50.9439570494;    // 51V
		case  439: return 51.9447730195;    // 52V
		case  440: return 52.944336734;     // 53V
		case  441: return 53.94643916;      // 54V
		case  442: return 54.9472410;       // 55V
		case  443: return 55.9504819;       // 56V
		case  444: return 56.9525224;       // 57V
		case  445: return 57.9567214;       // 58V
		case  446: return 58.9593917;       // 59V
		case  447: return 59.9643124;       // 60V
		case  448: return 60.9672596;       // 61V
		case  449: return 61.9726532;       // 62V
		case  450: return 62.9763943;       // 63V
		case  451: return 63.9826443;       // 64V
		case  452: return 64.9875054;       // 65V
		case  453: return 65.9939864;       // 66V
		case  454: return 42.0067043;       // 42Cr
		case  455: return 42.9975343;       // 43Cr
		case  456: return 43.9853632;       // 44Cr
		case  457: return 44.97905038;      // 45Cr
		case  458: return 45.96835921;      // 46Cr
		case  459: return 46.962897475;     // 47Cr
		case  460: return 47.954029179;     // 48Cr
		case  461: return 48.951333325;     // 49Cr
		case  462: return 49.9460418394;    // 50Cr
		case  463: return 50.9447650294;    // 51Cr
		case  464: return 51.9405062363;    // 52Cr
		case  465: return 52.9406481562;    // 53Cr
		case  466: return 53.9388791661;    // 54Cr
		case  467: return 54.9408384364;    // 55Cr
		case  468: return 55.940653120;     // 56Cr
		case  469: return 56.943613020;     // 57Cr
		case  470: return 57.9443522;       // 58Cr
		case  471: return 58.9485926;       // 59Cr
		case  472: return 59.9500823;       // 60Cr
		case  473: return 60.9544214;       // 61Cr
		case  474: return 61.9561016;       // 62Cr
		case  475: return 62.9616549;       // 63Cr
		case  476: return 63.9640832;       // 64Cr
		case  477: return 64.9699632;       // 65Cr
		case  478: return 65.9736654;       // 66Cr
		case  479: return 66.9801654;       // 67Cr
		case  480: return 67.9840375;       // 68Cr
		case  481: return 44.0071554;       // 44Mn
		case  482: return 44.9944943;       // 45Mn
		case  483: return 45.9860943;       // 46Mn
		case  484: return 46.97577534;      // 47Mn
		case  485: return 47.9685218;       // 48Mn
		case  486: return 48.95959511;      // 49Mn
		case  487: return 49.9542377895;    // 50Mn
		case  488: return 50.9482084794;    // 51Mn
		case  489: return 51.945563920;     // 52Mn
		case  490: return 52.9412888968;    // 53Mn
		case  491: return 53.940357612;     // 54Mn
		case  492: return 54.9380439148;    // 55Mn
		case  493: return 55.9389036949;    // 56Mn
		case  494: return 56.938286116;     // 57Mn
		case  495: return 57.940066629;     // 58Mn
		case  496: return 58.940391125;     // 59Mn
		case  497: return 59.943136625;     // 60Mn
		case  498: return 60.944452525;     // 61Mn
		case  499: return 61.9479516;       // 62Mn
		case  500: return 62.949664740;     // 63Mn
		case  501: return 63.953849438;     // 64Mn
		case  502: return 64.956019840;     // 65Mn
		case  503: return 65.96054712;      // 66Mn
		case  504: return 66.9642443;       // 67Mn
		case  505: return 67.9696254;       // 68Mn
		case  506: return 68.9736664;       // 69Mn
		case  507: return 69.9793775;       // 70Mn
		case  508: return 70.9836875;       // 71Mn
		case  509: return 45.0144243;       // 45Fe
		case  510: return 46.0006354;       // 46Fe
		case  511: return 46.9918554;       // 47Fe
		case  512: return 47.9802343;       // 48Fe
		case  513: return 48.97342926;      // 49Fe
		case  514: return 49.96297564;      // 50Fe
		case  515: return 50.956841096;     // 51Fe
		case  516: return 51.948113170;     // 52Fe
		case  517: return 52.945306418;     // 53Fe
		case  518: return 53.9396089953;    // 54Fe
		case  519: return 54.9382919951;    // 55Fe
		case  520: return 55.9349363349;    // 56Fe
		case  521: return 56.9353928449;    // 57Fe
		case  522: return 57.9332744353;    // 58Fe
		case  523: return 58.9348743454;    // 59Fe
		case  524: return 59.934071137;     // 60Fe
		case  525: return 60.936746228;     // 61Fe
		case  526: return 61.936791830;     // 62Fe
		case  527: return 62.940272746;     // 63Fe
		case  528: return 63.940987854;     // 64Fe
		case  529: return 64.945011573;     // 65Fe
		case  530: return 65.946250044;     // 66Fe
		case  531: return 66.9505423;       // 67Fe
		case  532: return 67.9529539;       // 68Fe
		case  533: return 68.9580743;       // 69Fe
		case  534: return 69.9610254;       // 70Fe
		case  535: return 70.9667264;       // 71Fe
		case  536: return 71.9698375;       // 72Fe
		case  537: return 72.9757275;       // 73Fe
		case  538: return 73.9793586;       // 74Fe
		case  539: return 47.0105786;       // 47Co
		case  540: return 48.0009386;       // 48Co
		case  541: return 48.9889175;       // 49Co
		case  542: return 49.9809164;       // 50Co
		case  543: return 50.97064752;      // 51Co
		case  544: return 51.9635121;       // 52Co
		case  545: return 52.954204119;     // 53Co
		case  546: return 53.9484598754;    // 54Co
		case  547: return 54.9419972057;    // 55Co
		case  548: return 55.9398388063;    // 56Co
		case  549: return 56.9362905766;    // 57Co
		case  550: return 57.935752113;     // 58Co
		case  551: return 58.9331942956;    // 59Co
		case  552: return 59.9338163056;    // 60Co
		case  553: return 60.9324766295;    // 61Co
		case  554: return 61.93405920;      // 62Co
		case  555: return 62.93360020;      // 63Co
		case  556: return 63.93581121;      // 64Co
		case  557: return 64.936462122;     // 65Co
		case  558: return 65.93944315;      // 66Co
		case  559: return 66.940609669;     // 67Co
		case  560: return 67.9442616;       // 68Co
		case  561: return 68.9461420;       // 69Co
		case  562: return 69.9496332;       // 70Co
		case  563: return 70.9523750;       // 71Co
		case  564: return 71.9572943;       // 72Co
		case  565: return 72.9603954;       // 73Co
		case  566: return 73.9651564;       // 74Co
		case  567: return 74.9687675;       // 75Co
		case  568: return 75.9741386;       // 76Co
		case  569: return 48.0176954;       // 48Ni
		case  570: return 49.0077086;       // 49Ni
		case  571: return 49.9947486;       // 50Ni
		case  572: return 50.9861186;       // 51Ni
		case  573: return 51.9748075;       // 52Ni
		case  574: return 52.96819027;      // 53Ni
		case  575: return 53.95789254;      // 54Ni
		case  576: return 54.9513306385;    // 55Ni
		case  577: return 55.9421285557;    // 56Ni
		case  578: return 56.9397921871;    // 57Ni
		case  579: return 57.9353424152;    // 58Ni
		case  580: return 58.9343462052;    // 59Ni
		case  581: return 59.9307858852;    // 60Ni
		case  582: return 60.9310555752;    // 61Ni
		case  583: return 61.9283453755;    // 62Ni
		case  584: return 62.9296696356;    // 63Ni
		case  585: return 63.9279668258;    // 64Ni
		case  586: return 64.9300851760;    // 65Ni
		case  587: return 65.929139315;     // 66Ni
		case  588: return 66.931569431;     // 67Ni
		case  589: return 67.931868832;     // 68Ni
		case  590: return 68.935610340;     // 69Ni
		case  591: return 69.936431323;     // 70Ni
		case  592: return 70.940519024;     // 71Ni
		case  593: return 71.941785924;     // 72Ni
		case  594: return 72.946206726;     // 73Ni
		case  595: return 73.9479843;       // 74Ni
		case  596: return 74.9525032;       // 75Ni
		case  597: return 75.9553354;       // 76Ni
		case  598: return 76.9605554;       // 77Ni
		case  599: return 77.9633686;       // 78Ni
		case  600: return 78.9702586;       // 79Ni
		case  601: return 51.9967186;       // 52Cu
		case  602: return 52.9845986;       // 53Cu
		case  603: return 53.9766654;       // 54Cu
		case  604: return 54.9660417;       // 55Cu
		case  605: return 55.9589521;       // 56Cu
		case  606: return 56.9492125066;    // 57Cu
		case  607: return 57.9445330570;    // 58Cu
		case  608: return 58.9394974867;    // 59Cu
		case  609: return 59.937364518;     // 60Cu
		case  610: return 60.933457610;     // 61Cu
		case  611: return 61.9325954175;    // 62Cu
		case  612: return 62.9295977256;    // 63Cu
		case  613: return 63.9297643456;    // 64Cu
		case  614: return 64.9277897071;    // 65Cu
		case  615: return 65.9288690372;    // 66Cu
		case  616: return 66.927730313;     // 67Cu
		case  617: return 67.929610917;     // 68Cu
		case  618: return 68.929429315;     // 69Cu
		case  619: return 69.932392112;     // 70Cu
		case  620: return 70.932676816;     // 71Cu
		case  621: return 71.935820315;     // 72Cu
		case  622: return 72.936674421;     // 73Cu
		case  623: return 73.939874966;     // 74Cu
		case  624: return 74.941522625;     // 75Cu
		case  625: return 75.945275072;     // 76Cu
		case  626: return 76.9479216;       // 77Cu
		case  627: return 77.9522354;       // 78Cu
		case  628: return 78.9550243;       // 79Cu
		case  629: return 79.9608964;       // 80Cu
		case  630: return 80.9658786;       // 81Cu
		case  631: return 81.9724486;       // 82Cu
		case  632: return 53.9920475;       // 54Zn
		case  633: return 54.9839875;       // 55Zn
		case  634: return 55.9725454;       // 56Zn
		case  635: return 56.9650622;       // 57Zn
		case  636: return 57.95459154;      // 58Zn
		case  637: return 58.9493126689;    // 59Zn
		case  638: return 59.9418421069;    // 60Zn
		case  639: return 60.93950717;      // 61Zn
		case  640: return 61.9343339773;    // 62Zn
		case  641: return 62.933211517;     // 63Zn
		case  642: return 63.9291420171;    // 64Zn
		case  643: return 64.9292407771;    // 65Zn
		case  644: return 65.9260338194;    // 66Zn
		case  645: return 66.9271277596;    // 67Zn
		case  646: return 67.9248445598;    // 68Zn
		case  647: return 68.926550710;     // 69Zn
		case  648: return 69.925319221;     // 70Zn
		case  649: return 70.927719628;     // 71Zn
		case  650: return 71.926842823;     // 72Zn
		case  651: return 72.929582620;     // 73Zn
		case  652: return 73.929407327;     // 74Zn
		case  653: return 74.932840221;     // 75Zn
		case  654: return 75.933115016;     // 76Zn
		case  655: return 76.936887221;     // 77Zn
		case  656: return 77.938289221;     // 78Zn
		case  657: return 78.942638124;     // 79Zn
		case  658: return 79.944552928;     // 80Zn
		case  659: return 80.950402654;     // 81Zn
		case  660: return 81.9542632;       // 82Zn
		case  661: return 82.9605654;       // 83Zn
		case  662: return 83.9652164;       // 84Zn
		case  663: return 84.9722675;       // 85Zn
		case  664: return 55.9953664;       // 56Ga
		case  665: return 56.9832032;       // 57Ga
		case  666: return 57.9747821;       // 58Ga
		case  667: return 58.9635318;       // 59Ga
		case  668: return 59.9572921;       // 60Ga
		case  669: return 60.94939941;      // 61Ga
		case  670: return 61.9441902575;    // 62Ga
		case  671: return 62.939294214;     // 63Ga
		case  672: return 63.936840415;     // 64Ga
		case  673: return 64.9327345988;    // 65Ga
		case  674: return 65.931589434;     // 66Ga
		case  675: return 66.928202513;     // 67Ga
		case  676: return 67.927980516;     // 68Ga
		case  677: return 68.925573513;     // 69Ga
		case  678: return 69.926021913;     // 70Ga
		case  679: return 70.9247025887;    // 71Ga
		case  680: return 71.9263674788;    // 72Ga
		case  681: return 72.925174718;     // 73Ga
		case  682: return 73.926945732;     // 74Ga
		case  683: return 74.926500226;     // 75Ga
		case  684: return 75.928827621;     // 76Ga
		case  685: return 76.929154326;     // 77Ga
		case  686: return 77.931608820;     // 78Ga
		case  687: return 78.932852320;     // 79Ga
		case  688: return 79.936420831;     // 80Ga
		case  689: return 80.938133835;     // 81Ga
		case  690: return 81.943176526;     // 82Ga
		case  691: return 82.947120328;     // 83Ga
		case  692: return 83.9524643;       // 84Ga
		case  693: return 84.9569932;       // 85Ga
		case  694: return 85.9630175;       // 86Ga
		case  695: return 86.9682486;       // 87Ga
		case  696: return 57.9917243;       // 58Ge
		case  697: return 58.9824932;       // 59Ge
		case  698: return 59.9703621;       // 60Ge
		case  699: return 60.9637932;       // 61Ge
		case  700: return 61.9550215;       // 62Ge
		case  701: return 62.94962840;      // 63Ge
		case  702: return 63.941689940;     // 64Ge
		case  703: return 64.939368123;     // 65Ge
		case  704: return 65.933862126;     // 66Ge
		case  705: return 66.932733950;     // 67Ge
		case  706: return 67.928095320;     // 68Ge
		case  707: return 68.927964514;     // 69Ge
		case  708: return 69.9242487590;    // 70Ge
		case  709: return 70.9249523390;    // 71Ge
		case  710: return 71.92207582681;   // 72Ge
		case  711: return 72.92345895661;   // 73Ge
		case  712: return 73.92117776113;   // 74Ge
		case  713: return 74.92285837055;   // 75Ge
		case  714: return 75.92140272619;   // 76Ge
		case  715: return 76.92354984357;   // 77Ge
		case  716: return 77.922852938;     // 78Ge
		case  717: return 78.92536040;      // 79Ge
		case  718: return 79.925350822;     // 80Ge
		case  719: return 80.928832922;     // 81Ge
		case  720: return 81.929774024;     // 82Ge
		case  721: return 82.934539126;     // 83Ge
		case  722: return 83.937575134;     // 84Ge
		case  723: return 84.942969740;     // 85Ge
		case  724: return 85.9465832;       // 86Ge
		case  725: return 86.9526843;       // 87Ge
		case  726: return 87.9569154;       // 88Ge
		case  727: return 88.9637964;       // 89Ge
		case  728: return 89.9686375;       // 90Ge
		case  729: return 59.9938843;       // 60As
		case  730: return 60.9811232;       // 61As
		case  731: return 61.9736132;       // 62As
		case  732: return 62.9639021;       // 63As
		case  733: return 63.9574333;       // 64As
		case  734: return 64.94961191;      // 65As
		case  735: return 65.944148861;     // 66As
		case  736: return 66.9392511148;    // 67As
		case  737: return 67.936774120;     // 68As
		case  738: return 68.93224634;      // 69As
		case  739: return 69.93092654;      // 70As
		case  740: return 70.927113845;     // 71As
		case  741: return 71.926752344;     // 72As
		case  742: return 72.923829141;     // 73As
		case  743: return 73.923928618;     // 74As
		case  744: return 74.9215945795;    // 75As
		case  745: return 75.9223920295;    // 76As
		case  746: return 76.920647618;     // 77As
		case  747: return 77.92182811;      // 78As
		case  748: return 78.920948458;     // 79As
		case  749: return 79.922474636;     // 80As
		case  750: return 80.922132329;     // 81As
		case  751: return 81.924741246;     // 82As
		case  752: return 82.925206930;     // 83As
		case  753: return 83.929303334;     // 84As
		case  754: return 84.932163733;     // 85As
		case  755: return 85.936701537;     // 86As
		case  756: return 86.940291732;     // 87As
		case  757: return 87.9455521;       // 88As
		case  758: return 88.9497632;       // 89As
		case  759: return 89.9556364;       // 90As
		case  760: return 90.9603964;       // 91As
		case  761: return 91.9667475;       // 92As
		case  762: return 63.9710954;       // 64Se
		case  763: return 64.9644064;       // 65Se
		case  764: return 65.9555932;       // 66Se
		case  765: return 66.94999472;      // 67Se
		case  766: return 67.9418252453;    // 68Se
		case  767: return 68.939414816;     // 69Se
		case  768: return 69.933515517;     // 70Se
		case  769: return 70.932209430;     // 71Se
		case  770: return 71.927140521;     // 72Se
		case  771: return 72.926754980;     // 73Se
		case  772: return 73.92247593415;   // 74Se
		case  773: return 74.92252287078;   // 75Se
		case  774: return 75.91921370417;   // 76Se
		case  775: return 76.91991415467;   // 77Se
		case  776: return 77.9173092820;    // 78Se
		case  777: return 78.9184992924;    // 79Se
		case  778: return 79.916521813;     // 80Se
		case  779: return 80.917993014;     // 81Se
		case  780: return 81.916699515;     // 82Se
		case  781: return 82.919118636;     // 83Se
		case  782: return 83.918466821;     // 84Se
		case  783: return 84.922260828;     // 85Se
		case  784: return 85.924311727;     // 86Se
		case  785: return 86.928688624;     // 87Se
		case  786: return 87.931417536;     // 88Se
		case  787: return 88.936669140;     // 89Se
		case  788: return 89.9401035;       // 90Se
		case  789: return 90.9459654;       // 91Se
		case  790: return 91.9498464;       // 92Se
		case  791: return 92.9562986;       // 93Se
		case  792: return 93.9604986;       // 94Se
		case  793: return 94.9673086;       // 95Se
		case  794: return 66.9646554;       // 67Br
		case  795: return 67.9587333;       // 68Br
		case  796: return 68.95049740;      // 69Br
		case  797: return 69.94479216;      // 70Br
		case  798: return 70.939342258;     // 71Br
		case  799: return 71.936588672;     // 72Br
		case  800: return 72.931671578;     // 73Br
		case  801: return 73.929910263;     // 74Br
		case  802: return 74.925810546;     // 75Br
		case  803: return 75.92454210;      // 76Br
		case  804: return 76.921379230;     // 77Br
		case  805: return 77.921145938;     // 78Br
		case  806: return 78.918337614;     // 79Br
		case  807: return 79.918529814;     // 80Br
		case  808: return 80.916289714;     // 81Br
		case  809: return 81.916803214;     // 82Br
		case  810: return 82.915175641;     // 83Br
		case  811: return 83.91649628;      // 84Br
		case  812: return 84.915645833;     // 85Br
		case  813: return 85.918805433;     // 86Br
		case  814: return 86.920674034;     // 87Br
		case  815: return 87.924083334;     // 88Br
		case  816: return 88.926704635;     // 89Br
		case  817: return 89.931292836;     // 90Br
		case  818: return 90.934398638;     // 91Br
		case  819: return 91.939631672;     // 92Br
		case  820: return 92.9431348;       // 93Br
		case  821: return 93.9489043;       // 94Br
		case  822: return 94.9530121;       // 95Br
		case  823: return 95.9590332;       // 96Br
		case  824: return 96.9634443;       // 97Br
		case  825: return 97.9694643;       // 98Br
		case  826: return 68.9651843;       // 69Kr
		case  827: return 69.9560421;       // 70Kr
		case  828: return 70.9502714;       // 71Kr
		case  829: return 71.942092486;     // 72Kr
		case  830: return 72.939289271;     // 73Kr
		case  831: return 73.933084022;     // 74Kr
		case  832: return 74.930945787;     // 75Kr
		case  833: return 75.925910343;     // 76Kr
		case  834: return 76.924670021;     // 77Kr
		case  835: return 77.9203649476;    // 78Kr
		case  836: return 78.920082938;     // 79Kr
		case  837: return 79.9163780875;    // 80Kr
		case  838: return 80.916591215;     // 81Kr
		case  839: return 81.9134827394;    // 82Kr
		case  840: return 82.9141271632;    // 83Kr
		case  841: return 83.911497728244;  // 84Kr
		case  842: return 84.912527321;     // 85Kr
		case  843: return 85.910610626941;  // 86Kr
		case  844: return 86.9133547626;    // 87Kr
		case  845: return 87.914447928;     // 88Kr
		case  846: return 88.917835523;     // 89Kr
		case  847: return 89.919527920;     // 90Kr
		case  848: return 90.923806324;     // 91Kr
		case  849: return 91.926173129;     // 92Kr
		case  850: return 92.931147227;     // 93Kr
		case  851: return 93.93414013;      // 94Kr
		case  852: return 94.93971120;      // 95Kr
		case  853: return 95.94301722;      // 96Kr
		case  854: return 96.9490914;       // 97Kr
		case  855: return 97.9524332;       // 98Kr
		case  856: return 98.9583954;       // 99Kr
		case  857: return 99.9623743;       // 100Kr
		case  858: return 100.9687354;      // 101Kr
		case  859: return 70.9653254;       // 71Rb
		case  860: return 71.9590854;       // 72Rb
		case  861: return 72.9505311;       // 73Rb
		case  862: return 73.944265932;     // 74Rb
		case  863: return 74.938573213;     // 75Rb
		case  864: return 75.935073010;     // 76Rb
		case  865: return 76.930401614;     // 77Rb
		case  866: return 77.928141935;     // 78Rb
		case  867: return 78.923989923;     // 79Rb
		case  868: return 79.922516420;     // 80Rb
		case  869: return 80.918993953;     // 81Rb
		case  870: return 81.918209032;     // 82Rb
		case  871: return 82.915114225;     // 83Rb
		case  872: return 83.914375224;     // 84Rb
		case  873: return 84.911789737954;  // 85Rb
		case  874: return 85.9111674321;    // 86Rb
		case  875: return 86.909180531060;  // 87Rb
		case  876: return 87.9113155917;    // 88Rb
		case  877: return 88.912278359;     // 89Rb
		case  878: return 89.914798570;     // 90Rb
		case  879: return 90.916537284;     // 91Rb
		case  880: return 91.919728466;     // 92Rb
		case  881: return 92.922039384;     // 93Rb
		case  882: return 93.926394822;     // 94Rb
		case  883: return 94.92926022;      // 95Rb
		case  884: return 95.934133436;     // 96Rb
		case  885: return 96.937177121;     // 97Rb
		case  886: return 97.941686937;     // 98Rb
		case  887: return 98.9450312;       // 99Rb
		case  888: return 99.9500321;       // 100Rb
		case  889: return 100.9540423;      // 101Rb
		case  890: return 101.9595232;      // 102Rb
		case  891: return 102.9639243;      // 103Rb
		case  892: return 72.9657043;       // 73Sr
		case  893: return 73.9561711;       // 74Sr
		case  894: return 74.9499524;       // 75Sr
		case  895: return 75.94176337;      // 76Sr
		case  896: return 76.937945585;     // 77Sr
		case  897: return 77.932180080;     // 78Sr
		case  898: return 78.929707790;     // 79Sr
		case  899: return 79.924517537;     // 80Sr
		case  900: return 80.923211434;     // 81Sr
		case  901: return 81.918399964;     // 82Sr
		case  902: return 82.917554473;     // 83Sr
		case  903: return 83.913419113;     // 84Sr
		case  904: return 84.912932030;     // 85Sr
		case  905: return 85.909260612;     // 86Sr
		case  906: return 86.908877512;     // 87Sr
		case  907: return 87.905612512;     // 88Sr
		case  908: return 88.907451112;     // 89Sr
		case  909: return 89.907730028;     // 90Sr
		case  910: return 90.910195461;     // 91Sr
		case  911: return 91.911038237;     // 92Sr
		case  912: return 92.914024281;     // 93Sr
		case  913: return 93.915355618;     // 94Sr
		case  914: return 94.919352963;     // 95Sr
		case  915: return 95.921706693;     // 96Sr
		case  916: return 96.926374036;     // 97Sr
		case  917: return 97.928688840;     // 98Sr
		case  918: return 98.932890738;     // 99Sr
		case  919: return 99.93577010;      // 100Sr
		case  920: return 100.94035286;     // 101Sr
		case  921: return 101.94379175;     // 102Sr
		case  922: return 102.9490921;      // 103Sr
		case  923: return 103.9526532;      // 104Sr
		case  924: return 104.9585554;      // 105Sr
		case  925: return 105.9626564;      // 106Sr
		case  926: return 106.9689775;      // 107Sr
		case  927: return 75.9585654;       // 76Y
		case  928: return 76.94978165;      // 77Y
		case  929: return 77.9436143;       // 78Y
		case  930: return 78.9373548;       // 79Y
		case  931: return 79.934356167;     // 80Y
		case  932: return 80.929455658;     // 81Y
		case  933: return 81.926931459;     // 82Y
		case  934: return 82.92248520;      // 83Y
		case  935: return 83.920672146;     // 84Y
		case  936: return 84.91643320;      // 85Y
		case  937: return 85.91488615;      // 86Y
		case  938: return 86.910876117;     // 87Y
		case  939: return 87.909501620;     // 88Y
		case  940: return 88.905840324;     // 89Y
		case  941: return 89.907143924;     // 90Y
		case  942: return 90.907297428;     // 91Y
		case  943: return 91.908945199;     // 92Y
		case  944: return 92.90957811;      // 93Y
		case  945: return 93.911590669;     // 94Y
		case  946: return 94.912816174;     // 95Y
		case  947: return 95.915896869;     // 96Y
		case  948: return 96.918274175;     // 97Y
		case  949: return 97.922382188;     // 98Y
		case  950: return 98.924148074;     // 99Y
		case  951: return 99.92771512;      // 100Y
		case  952: return 100.930147779;    // 101Y
		case  953: return 101.934327744;    // 102Y
		case  954: return 102.93724312;     // 103Y
		case  955: return 103.9419643;      // 104Y
		case  956: return 104.9454454;      // 105Y
		case  957: return 105.9505654;      // 106Y
		case  958: return 106.9545254;      // 107Y
		case  959: return 107.9599664;      // 108Y
		case  960: return 108.9643675;      // 109Y
		case  961: return 77.9556654;       // 78Zr
		case  962: return 78.9494843;       // 79Zr
		case  963: return 79.940416;        // 80Zr
		case  964: return 80.9373118;       // 81Zr
		case  965: return 81.9313522;       // 82Zr
		case  966: return 82.929242169;     // 83Zr
		case  967: return 83.923326959;     // 84Zr
		case  968: return 84.921444469;     // 85Zr
		case  969: return 85.916297238;     // 86Zr
		case  970: return 86.914818045;     // 87Zr
		case  971: return 87.910221358;     // 88Zr
		case  972: return 88.908881437;     // 89Zr
		case  973: return 89.904697720;     // 90Zr
		case  974: return 90.905639620;     // 91Zr
		case  975: return 91.905034720;     // 92Zr
		case  976: return 92.906469920;     // 93Zr
		case  977: return 93.906310820;     // 94Zr
		case  978: return 94.908038519;     // 95Zr
		case  979: return 95.908271421;     // 96Zr
		case  980: return 96.910951221;     // 97Zr
		case  981: return 97.912728993;     // 98Zr
		case  982: return 98.91666711;      // 99Zr
		case  983: return 99.918000689;     // 100Zr
		case  984: return 100.921448091;    // 101Zr
		case  985: return 101.923140997;    // 102Zr
		case  986: return 102.92719110;     // 103Zr
		case  987: return 103.92943610;     // 104Zr
		case  988: return 104.93400813;     // 105Zr
		case  989: return 105.9367621;      // 106Zr
		case  990: return 106.9417432;      // 107Zr
		case  991: return 107.9448743;      // 108Zr
		case  992: return 108.9504154;      // 109Zr
		case  993: return 109.9539664;      // 110Zr
		case  994: return 110.9596875;      // 111Zr
		case  995: return 111.9637075;      // 112Zr
		case  996: return 80.9496043;       // 81Nb
		case  997: return 81.9439632;       // 82Nb
		case  998: return 82.9372932;       // 83Nb
		case  999: return 83.9344932;       // 84Nb
		case 1000: return 84.928845844;     // 85Nb
		case 1001: return 85.925782859;     // 86Nb
		case 1002: return 86.920693773;     // 87Nb
		case 1003: return 87.91822261;      // 88Nb
		case 1004: return 88.91344525;      // 89Nb
		case 1005: return 89.911258438;     // 90Nb
		case 1006: return 90.906989737;     // 91Nb
		case 1007: return 91.907188126;     // 92Nb
		case 1008: return 92.906373020;     // 93Nb
		case 1009: return 93.907278820;     // 94Nb
		case 1010: return 94.9068324071;    // 95Nb
		case 1011: return 95.908097335;     // 96Nb
		case 1012: return 96.908095919;     // 97Nb
		case 1013: return 97.910326558;     // 98Nb
		case 1014: return 98.91161313;      // 99Nb
		case 1015: return 99.914327688;     // 100Nb
		case 1016: return 100.915310342;    // 101Nb
		case 1017: return 101.918077235;    // 102Nb
		case 1018: return 102.919457244;    // 103Nb
		case 1019: return 103.922892537;    // 104Nb
		case 1020: return 104.924946545;    // 105Nb
		case 1021: return 105.928931746;    // 106Nb
		case 1022: return 106.931593787;    // 107Nb
		case 1023: return 107.936074888;    // 108Nb
		case 1024: return 108.9392256;      // 109Nb
		case 1025: return 109.9440321;      // 110Nb
		case 1026: return 110.9475332;      // 111Nb
		case 1027: return 111.9524732;      // 112Nb
		case 1028: return 112.9565143;      // 113Nb
		case 1029: return 113.9620154;      // 114Nb
		case 1030: return 114.9663454;      // 115Nb
		case 1031: return 82.9498843;       // 83Mo
		case 1032: return 83.9414943;       // 84Mo
		case 1033: return 84.93826117;      // 85Mo
		case 1034: return 85.931174840;     // 86Mo
		case 1035: return 86.928196231;     // 87Mo
		case 1036: return 87.921967841;     // 88Mo
		case 1037: return 88.919468242;     // 89Mo
		case 1038: return 89.913930938;     // 90Mo
		case 1039: return 90.911745367;     // 91Mo
		case 1040: return 91.9068079684;    // 92Mo
		case 1041: return 92.9068095884;    // 93Mo
		case 1042: return 93.9050849048;    // 94Mo
		case 1043: return 94.9058387747;    // 95Mo
		case 1044: return 95.9046761247;    // 96Mo
		case 1045: return 96.9060181249;    // 97Mo
		case 1046: return 97.9054048249;    // 98Mo
		case 1047: return 98.9077085152;    // 99Mo
		case 1048: return 99.907471811;     // 100Mo
		case 1049: return 100.910341411;    // 101Mo
		case 1050: return 101.910283491;    // 102Mo
		case 1051: return 102.91307910;     // 103Mo
		case 1052: return 103.913734498;    // 104Mo
		case 1053: return 104.91696910;     // 105Mo
		case 1054: return 105.91825910;     // 106Mo
		case 1055: return 106.92210610;     // 107Mo
		case 1056: return 107.92403310;     // 108Mo
		case 1057: return 108.92842412;     // 109Mo
		case 1058: return 109.93070426;     // 110Mo
		case 1059: return 110.93565414;     // 111Mo
		case 1060: return 111.9383121;      // 112Mo
		case 1061: return 112.9433532;      // 113Mo
		case 1062: return 113.9465332;      // 114Mo
		case 1063: return 114.9519643;      // 115Mo
		case 1064: return 115.9554554;      // 116Mo
		case 1065: return 116.9611754;      // 117Mo
		case 1066: return 84.9505843;       // 85Tc
		case 1067: return 85.9449332;       // 86Tc
		case 1068: return 86.938067245;     // 87Tc
		case 1069: return 87.9337816;       // 88Tc
		case 1070: return 88.927648741;     // 89Tc
		case 1071: return 89.924073911;     // 90Tc
		case 1072: return 90.918425425;     // 91Tc
		case 1073: return 91.915269833;     // 92Tc
		case 1074: return 92.910246014;     // 93Tc
		case 1075: return 93.909653644;     // 94Tc
		case 1076: return 94.907653655;     // 95Tc
		case 1077: return 95.907868055;     // 96Tc
		case 1078: return 96.906366740;     // 97Tc
		case 1079: return 97.907212436;     // 98Tc
		case 1080: return 98.906250810;     // 99Tc
		case 1081: return 99.907653915;     // 100Tc
		case 1082: return 100.90730926;     // 101Tc
		case 1083: return 101.909209799;    // 102Tc
		case 1084: return 102.90917611;     // 103Tc
		case 1085: return 103.91142527;     // 104Tc
		case 1086: return 104.91165538;     // 105Tc
		case 1087: return 105.91435813;     // 106Tc
		case 1088: return 106.915460693;    // 107Tc
		case 1089: return 107.918495794;    // 108Tc
		case 1090: return 108.92025610;     // 109Tc
		case 1091: return 109.92374410;     // 110Tc
		case 1092: return 110.92590111;     // 111Tc
		case 1093: return 111.929945860;    // 112Tc
		case 1094: return 112.932569036;    // 113Tc
		case 1095: return 113.9369111;      // 114Tc
		case 1096: return 114.9399821;      // 115Tc
		case 1097: return 115.9447632;      // 116Tc
		case 1098: return 116.9480643;      // 117Tc
		case 1099: return 117.9529943;      // 118Tc
		case 1100: return 118.9566654;      // 119Tc
		case 1101: return 119.9618754;      // 120Tc
		case 1102: return 86.9506943;       // 87Ru
		case 1103: return 87.9416032;       // 88Ru
		case 1104: return 88.9376232;       // 89Ru
		case 1105: return 89.930344440;     // 90Ru
		case 1106: return 90.926741924;     // 91Ru
		case 1107: return 91.920234429;     // 92Ru
		case 1108: return 92.917104422;     // 93Ru
		case 1109: return 93.911342934;     // 94Ru
		case 1110: return 94.91040610;      // 95Ru
		case 1111: return 95.9075902549;    // 96Ru
		case 1112: return 96.907547130;     // 97Ru
		case 1113: return 97.905286869;     // 98Ru
		case 1114: return 98.905934111;     // 99Ru
		case 1115: return 99.904214311;     // 100Ru
		case 1116: return 100.905576912;    // 101Ru
		case 1117: return 101.904344112;    // 102Ru
		case 1118: return 102.906318612;    // 103Ru
		case 1119: return 103.905427528;    // 104Ru
		case 1120: return 104.907747628;    // 105Ru
		case 1121: return 105.907329158;    // 106Ru
		case 1122: return 106.909972093;    // 107Ru
		case 1123: return 107.910188093;    // 108Ru
		case 1124: return 108.913326096;    // 109Ru
		case 1125: return 109.914040796;    // 110Ru
		case 1126: return 110.91757010;     // 111Ru
		case 1127: return 111.91880910;     // 112Ru
		case 1128: return 112.92284439;     // 113Ru
		case 1129: return 113.924613638;    // 114Ru
		case 1130: return 114.92882071;     // 115Ru
		case 1131: return 115.931219240;    // 116Ru
		case 1132: return 116.9361063;      // 117Ru
		case 1133: return 117.9385332;      // 118Ru
		case 1134: return 118.9435732;      // 119Ru
		case 1135: return 119.9463143;      // 120Ru
		case 1136: return 120.9516443;      // 121Ru
		case 1137: return 121.9544754;      // 122Ru
		case 1138: return 122.9598954;      // 123Ru
		case 1139: return 123.9630564;      // 124Ru
		case 1140: return 88.9505839;       // 89Rh
		case 1141: return 89.9442243;       // 90Rh
		case 1142: return 90.9368843;       // 91Rh
		case 1143: return 91.932367747;     // 92Rh
		case 1144: return 92.925912828;     // 93Rh
		case 1145: return 93.921730536;     // 94Rh
		case 1146: return 94.915897942;     // 95Rh
		case 1147: return 95.91445311;      // 96Rh
		case 1148: return 96.91132938;      // 97Rh
		case 1149: return 97.91070813;      // 98Rh
		case 1150: return 98.908128273;     // 99Rh
		case 1151: return 99.90811719;      // 100Rh
		case 1152: return 100.906160663;    // 101Rh
		case 1153: return 101.906837450;    // 102Rh
		case 1154: return 102.905498026;    // 103Rh
		case 1155: return 103.906649226;    // 104Rh
		case 1156: return 104.905688527;    // 105Rh
		case 1157: return 105.907286858;    // 106Rh
		case 1158: return 106.90674813;     // 107Rh
		case 1159: return 107.90871415;     // 108Rh
		case 1160: return 108.908748843;    // 109Rh
		case 1161: return 109.91107919;     // 110Rh
		case 1162: return 110.911642374;    // 111Rh
		case 1163: return 111.91440347;     // 112Rh
		case 1164: return 112.915439377;    // 113Rh
		case 1165: return 113.91871877;     // 114Rh
		case 1166: return 114.920311678;    // 115Rh
		case 1167: return 115.92405976;     // 116Rh
		case 1168: return 116.926035495;    // 117Rh
		case 1169: return 117.93034026;     // 118Rh
		case 1170: return 118.93255710;     // 119Rh
		case 1171: return 119.9368621;      // 120Rh
		case 1172: return 120.9394232;      // 121Rh
		case 1173: return 121.9439932;      // 122Rh
		case 1174: return 122.9468543;      // 123Rh
		case 1175: return 123.9515143;      // 124Rh
		case 1176: return 124.9546954;      // 125Rh
		case 1177: return 125.9594654;      // 126Rh
		case 1178: return 90.9503254;       // 91Pd
		case 1179: return 91.9408854;       // 92Pd
		case 1180: return 92.9365143;       // 93Pd
		case 1181: return 93.929037646;     // 94Pd
		case 1182: return 94.924889833;     // 95Pd
		case 1183: return 95.918215145;     // 96Pd
		case 1184: return 96.916472052;     // 97Pd
		case 1185: return 97.912698351;     // 98Pd
		case 1186: return 98.911774854;     // 99Pd
		case 1187: return 99.90850519;      // 100Pd
		case 1188: return 100.908286449;    // 101Pd
		case 1189: return 101.905602228;    // 102Pd
		case 1190: return 102.906080927;    // 103Pd
		case 1191: return 103.904030514;    // 104Pd
		case 1192: return 104.905079612;    // 105Pd
		case 1193: return 105.903480412;    // 106Pd
		case 1194: return 106.905128213;    // 107Pd
		case 1195: return 107.903891612;    // 108Pd
		case 1196: return 108.905950412;    // 109Pd
		case 1197: return 109.9051722075;   // 110Pd
		case 1198: return 110.9076896886;   // 111Pd
		case 1199: return 111.907329770;    // 112Pd
		case 1200: return 112.910261075;    // 113Pd
		case 1201: return 113.910368675;    // 114Pd
		case 1202: return 114.91365915;     // 115Pd
		case 1203: return 115.914297077;    // 116Pd
		case 1204: return 116.917954778;    // 117Pd
		case 1205: return 117.919066727;    // 118Pd
		case 1206: return 118.923340289;    // 119Pd
		case 1207: return 119.924551125;    // 120Pd
		case 1208: return 120.928950336;    // 121Pd
		case 1209: return 121.93063221;     // 122Pd
		case 1210: return 122.9351421;      // 123Pd
		case 1211: return 123.9371432;      // 124Pd
		case 1212: return 124.9417943;      // 125Pd
		case 1213: return 125.9441654;      // 126Pd
		case 1214: return 126.9490754;      // 127Pd
		case 1215: return 127.9518364;      // 128Pd
		case 1216: return 92.9503354;       // 93Ag
		case 1217: return 93.9437369;       // 94Ag
		case 1218: return 94.9360243;       // 95Ag
		case 1219: return 95.93074497;      // 96Ag
		case 1220: return 96.9239712;       // 97Ag
		case 1221: return 97.92156035;      // 98Ag
		case 1222: return 98.917645867;     // 99Ag
		case 1223: return 99.916115454;     // 100Ag
		case 1224: return 100.912684052;    // 101Ag
		case 1225: return 101.911704788;    // 102Ag
		case 1226: return 102.908963141;    // 103Ag
		case 1227: return 103.908623945;    // 104Ag
		case 1228: return 104.906525649;    // 105Ag
		case 1229: return 105.906663632;    // 106Ag
		case 1230: return 106.905091626;    // 107Ag
		case 1231: return 107.905950326;    // 108Ag
		case 1232: return 108.904755314;    // 109Ag
		case 1233: return 109.906110214;    // 110Ag
		case 1234: return 110.905295916;    // 111Ag
		case 1235: return 111.907048626;    // 112Ag
		case 1236: return 112.90657318;     // 113Ag
		case 1237: return 113.908823049;    // 114Ag
		case 1238: return 114.90876720;     // 115Ag
		case 1239: return 115.911386835;    // 116Ag
		case 1240: return 116.91177415;     // 117Ag
		case 1241: return 117.914595527;    // 118Ag
		case 1242: return 118.91557016;     // 119Ag
		case 1243: return 119.918784848;    // 120Ag
		case 1244: return 120.92012513;     // 121Ag
		case 1245: return 121.92366441;     // 122Ag
		case 1246: return 122.92533733;     // 123Ag
		case 1247: return 123.9289327;      // 124Ag
		case 1248: return 124.9310564;      // 125Ag
		case 1249: return 125.9347521;      // 126Ag
		case 1250: return 126.9371121;      // 127Ag
		case 1251: return 127.9410632;      // 128Ag
		case 1252: return 128.9439532;      // 129Ag
		case 1253: return 129.9507036;      // 130Ag
		case 1254: return 94.9499454;       // 95Cd
		case 1255: return 95.9403443;       // 96Cd
		case 1256: return 96.9351032;       // 97Cd
		case 1257: return 97.92738956;      // 98Cd
		case 1258: return 98.924925817;     // 99Cd
		case 1259: return 99.920348818;     // 100Cd
		case 1260: return 100.918586216;    // 101Cd
		case 1261: return 101.914482018;    // 102Cd
		case 1262: return 102.913416519;    // 103Cd
		case 1263: return 103.909856418;    // 104Cd
		case 1264: return 104.909463915;    // 105Cd
		case 1265: return 105.906459912;    // 106Cd
		case 1266: return 106.906612118;    // 107Cd
		case 1267: return 107.904183412;    // 108Cd
		case 1268: return 108.904986717;    // 109Cd
		case 1269: return 109.9030066161;   // 110Cd
		case 1270: return 110.9041828761;   // 111Cd
		case 1271: return 111.9027628760;   // 112Cd
		case 1272: return 112.9044081345;   // 113Cd
		case 1273: return 113.9033650943;   // 114Cd
		case 1274: return 114.9054375177;   // 115Cd
		case 1275: return 115.9047631517;   // 116Cd
		case 1276: return 116.907226011;    // 117Cd
		case 1277: return 117.90692221;     // 118Cd
		case 1278: return 118.90984740;     // 119Cd
		case 1279: return 119.909868140;    // 120Cd
		case 1280: return 120.912963721;    // 121Cd
		case 1281: return 121.913459125;    // 122Cd
		case 1282: return 122.916892529;    // 123Cd
		case 1283: return 123.917657432;    // 124Cd
		case 1284: return 124.921257631;    // 125Cd
		case 1285: return 125.922429127;    // 126Cd
		case 1286: return 126.92647214;     // 127Cd
		case 1287: return 127.927812978;    // 128Cd
		case 1288: return 128.9318221;      // 129Cd
		case 1289: return 129.9339418;      // 130Cd
		case 1290: return 130.9406021;      // 131Cd
		case 1291: return 131.9460421;      // 132Cd
		case 1292: return 132.9528532;      // 133Cd
		case 1293: return 96.9493454;       // 97In
		case 1294: return 97.9421421;       // 98In
		case 1295: return 98.9341121;       // 99In
		case 1296: return 99.9309620;       // 100In
		case 1297: return 100.9263432;      // 101In
		case 1298: return 101.924107149;    // 102In
		case 1299: return 102.919881998;    // 103In
		case 1300: return 103.918214562;    // 104In
		case 1301: return 104.91450211;     // 105In
		case 1302: return 105.91346413;     // 106In
		case 1303: return 106.91029012;     // 107In
		case 1304: return 107.909693593;    // 108In
		case 1305: return 108.907151443;    // 109In
		case 1306: return 109.90717012;     // 110In
		case 1307: return 110.905108538;    // 111In
		case 1308: return 111.905537746;    // 112In
		case 1309: return 112.9040618491;   // 113In
		case 1310: return 113.9049179194;   // 114In
		case 1311: return 114.90387877612;  // 115In
		case 1312: return 115.9052599924;   // 116In
		case 1313: return 116.904515752;    // 117In
		case 1314: return 117.906356683;    // 118In
		case 1315: return 118.905850778;    // 119In
		case 1316: return 119.90796743;     // 120In
		case 1317: return 120.90785129;     // 121In
		case 1318: return 121.91028154;     // 122In
		case 1319: return 122.91043421;     // 123In
		case 1320: return 123.91318233;     // 124In
		case 1321: return 124.91360529;     // 125In
		case 1322: return 125.91650729;     // 126In
		case 1323: return 126.91744623;     // 127In
		case 1324: return 127.9204016;      // 128In
		case 1325: return 128.921805329;    // 129In
		case 1326: return 129.92497741;     // 130In
		case 1327: return 130.926971529;    // 131In
		case 1328: return 131.93300164;     // 132In
		case 1329: return 132.9383121;      // 133In
		case 1330: return 133.9445432;      // 134In
		case 1331: return 134.9500543;      // 135In
		case 1332: return 98.9485354;       // 99Sn
		case 1333: return 99.9385032;       // 100Sn
		case 1334: return 100.9352632;      // 101Sn
		case 1335: return 101.9302911;      // 102Sn
		case 1336: return 102.92810576;     // 103Sn
		case 1337: return 103.923105262;    // 104Sn
		case 1338: return 104.921268443;    // 105Sn
		case 1339: return 105.916957455;    // 106Sn
		case 1340: return 106.915713757;    // 107Sn
		case 1341: return 107.911894358;    // 108Sn
		case 1342: return 108.911292185;    // 109Sn
		case 1343: return 109.90784515;     // 110Sn
		case 1344: return 110.907740158;    // 111Sn
		case 1345: return 111.9048238761;   // 112Sn
		case 1346: return 112.905175718;    // 113Sn
		case 1347: return 113.902782710;    // 114Sn
		case 1348: return 114.90334469916;  // 115Sn
		case 1349: return 115.9017428010;   // 116Sn
		case 1350: return 116.9029539852;   // 117Sn
		case 1351: return 117.9016065754;   // 118Sn
		case 1352: return 118.9033111778;   // 119Sn
		case 1353: return 119.9022016397;   // 120Sn
		case 1354: return 120.904242610;    // 121Sn
		case 1355: return 121.903443826;    // 122Sn
		case 1356: return 122.905725226;    // 123Sn
		case 1357: return 123.905276611;    // 124Sn
		case 1358: return 124.907786411;    // 125Sn
		case 1359: return 125.90765911;     // 126Sn
		case 1360: return 126.91039011;     // 127Sn
		case 1361: return 127.91050719;     // 128Sn
		case 1362: return 128.91346521;     // 129Sn
		case 1363: return 129.913973823;    // 130Sn
		case 1364: return 130.917045065;    // 131Sn
		case 1365: return 131.917826731;    // 132Sn
		case 1366: return 132.923913426;    // 133Sn
		case 1367: return 133.928682135;    // 134Sn
		case 1368: return 134.934908633;    // 135Sn
		case 1369: return 135.9399943;      // 136Sn
		case 1370: return 136.9465554;      // 137Sn
		case 1371: return 137.9518464;      // 138Sn
		case 1372: return 102.9396932;      // 103Sb
		case 1373: return 103.9364813;      // 104Sb
		case 1374: return 104.93127623;     // 105Sb
		case 1375: return 105.928638080;    // 106Sb
		case 1376: return 106.924150645;    // 107Sb
		case 1377: return 107.922226759;    // 108Sb
		case 1378: return 108.918141157;    // 109Sb
		case 1379: return 109.916854364;    // 110Sb
		case 1380: return 110.913218295;    // 111Sb
		case 1381: return 111.91240019;     // 112Sb
		case 1382: return 112.90937518;     // 113Sb
		case 1383: return 113.90929023;     // 114Sb
		case 1384: return 114.90659817;     // 115Sb
		case 1385: return 115.906793155;    // 116Sb
		case 1386: return 116.904841591;    // 117Sb
		case 1387: return 117.905532132;    // 118Sb
		case 1388: return 118.903945583;    // 119Sb
		case 1389: return 119.905079477;    // 120Sb
		case 1390: return 120.903812030;    // 121Sb
		case 1391: return 121.905169930;    // 122Sb
		case 1392: return 122.904213223;    // 123Sb
		case 1393: return 123.905935023;    // 124Sb
		case 1394: return 124.905253028;    // 125Sb
		case 1395: return 125.90725334;     // 126Sb
		case 1396: return 126.906924355;    // 127Sb
		case 1397: return 127.90914621;     // 128Sb
		case 1398: return 128.90914723;     // 129Sb
		case 1399: return 129.91166215;     // 130Sb
		case 1400: return 130.911988823;    // 131Sb
		case 1401: return 131.914507729;    // 132Sb
		case 1402: return 132.915273234;    // 133Sb
		case 1403: return 133.920535718;    // 134Sb
		case 1404: return 134.925185131;    // 135Sb
		case 1405: return 135.930745968;    // 136Sb
		case 1406: return 136.9355532;      // 137Sb
		case 1407: return 137.9414532;      // 138Sb
		case 1408: return 138.9465543;      // 139Sb
		case 1409: return 139.9528364;      // 140Sb
		case 1410: return 104.9433032;      // 105Te
		case 1411: return 105.9375011;      // 106Te
		case 1412: return 106.93501276;     // 107Te
		case 1413: return 107.929380558;    // 108Te
		case 1414: return 108.927304547;    // 109Te
		case 1415: return 109.922458171;    // 110Te
		case 1416: return 110.921000669;    // 111Te
		case 1417: return 111.916727990;    // 112Te
		case 1418: return 112.91589130;     // 113Te
		case 1419: return 113.91208930;     // 114Te
		case 1420: return 114.91190230;     // 115Te
		case 1421: return 115.90846030;     // 116Te
		case 1422: return 116.90864614;     // 117Te
		case 1423: return 117.90585420;     // 118Te
		case 1424: return 118.906407185;    // 119Te
		case 1425: return 119.904059333;    // 120Te
		case 1426: return 120.90494428;     // 121Te
		case 1427: return 121.903043516;    // 122Te
		case 1428: return 122.904269816;    // 123Te
		case 1429: return 123.902817116;    // 124Te
		case 1430: return 124.904429916;    // 125Te
		case 1431: return 125.903310916;    // 126Te
		case 1432: return 126.905225716;    // 127Te
		case 1433: return 127.9044612893;   // 128Te
		case 1434: return 128.9065964693;   // 129Te
		case 1435: return 129.90622274812;  // 130Te
		case 1436: return 130.90852221365;  // 131Te
		case 1437: return 131.908546737;    // 132Te
		case 1438: return 132.910968839;    // 133Te
		case 1439: return 133.911394030;    // 134Te
		case 1440: return 134.916555729;    // 135Te
		case 1441: return 135.920100626;    // 136Te
		case 1442: return 136.925598927;    // 137Te
		case 1443: return 137.929472247;    // 138Te
		case 1444: return 138.935367238;    // 139Te
		case 1445: return 139.93949930;     // 140Te
		case 1446: return 140.9458043;      // 141Te
		case 1447: return 141.9502254;      // 142Te
		case 1448: return 142.9567654;      // 143Te
		case 1449: return 106.9467832;      // 107I
		case 1450: return 107.9434814;      // 108I
		case 1451: return 108.938085361;    // 109I
		case 1452: return 109.93508954;     // 110I
		case 1453: return 110.930269251;    // 111I
		case 1454: return 111.92800511;     // 112I
		case 1455: return 112.923650186;    // 113I
		case 1456: return 113.9218532;      // 114I
		case 1457: return 114.91804831;     // 115I
		case 1458: return 115.9168110;      // 116I
		case 1459: return 116.91364828;     // 117I
		case 1460: return 117.91307421;     // 118I
		case 1461: return 118.91007430;     // 119I
		case 1462: return 119.91008716;     // 120I
		case 1463: return 120.907405158;    // 121I
		case 1464: return 121.907588856;    // 122I
		case 1465: return 122.905588540;    // 123I
		case 1466: return 123.906209026;    // 124I
		case 1467: return 124.904629416;    // 125I
		case 1468: return 125.905623341;    // 126I
		case 1469: return 126.904471939;    // 127I
		case 1470: return 127.905808639;    // 128I
		case 1471: return 128.904983734;    // 129I
		case 1472: return 129.906670234;    // 130I
		case 1473: return 130.9061263069;   // 131I
		case 1474: return 131.907993544;    // 132I
		case 1475: return 132.907797050;    // 133I
		case 1476: return 133.909758859;    // 134I
		case 1477: return 134.910048858;    // 135I
		case 1478: return 135.91460415;     // 136I
		case 1479: return 136.918028290;    // 137I
		case 1480: return 137.922726464;    // 138I
		case 1481: return 138.92650631;     // 139I
		case 1482: return 139.9317320;      // 140I
		case 1483: return 140.9356921;      // 141I
		case 1484: return 141.9412040;      // 142I
		case 1485: return 142.9456532;      // 143I
		case 1486: return 143.9513943;      // 144I
		case 1487: return 144.9560554;      // 145I
		case 1488: return 108.9504332;      // 109Xe
		case 1489: return 109.9442611;      // 110Xe
		case 1490: return 110.94160793;     // 111Xe
		case 1491: return 111.935559089;    // 112Xe
		case 1492: return 112.933221773;    // 113Xe
		case 1493: return 113.92798012;     // 114Xe
		case 1494: return 114.92629413;     // 115Xe
		case 1495: return 115.92158114;     // 116Xe
		case 1496: return 116.92035911;     // 117Xe
		case 1497: return 117.91617911;     // 118Xe
		case 1498: return 118.91541111;     // 119Xe
		case 1499: return 119.91178413;     // 120Xe
		case 1500: return 120.91145311;     // 121Xe
		case 1501: return 121.90836812;     // 122Xe
		case 1502: return 122.90848210;     // 123Xe
		case 1503: return 123.905892019;    // 124Xe
		case 1504: return 124.906394420;    // 125Xe
		case 1505: return 125.904298338;    // 126Xe
		case 1506: return 126.905182944;    // 127Xe
		case 1507: return 127.903531011;    // 128Xe
		case 1508: return 128.904780861160; // 129Xe
		case 1509: return 129.90350934910;  // 130Xe
		case 1510: return 130.9050840624;   // 131Xe
		case 1511: return 131.904155085656; // 132Xe
		case 1512: return 132.905910826;    // 133Xe
		case 1513: return 133.9053946690;   // 134Xe
		case 1514: return 134.907227845;    // 135Xe
		case 1515: return 135.90721448411;  // 136Xe
		case 1516: return 136.9115577811;   // 137Xe
		case 1517: return 137.914146330;    // 138Xe
		case 1518: return 138.918792223;    // 139Xe
		case 1519: return 139.921645825;    // 140Xe
		case 1520: return 140.926787231;    // 141Xe
		case 1521: return 141.929973129;    // 142Xe
		case 1522: return 142.935369650;    // 143Xe
		case 1523: return 143.938945157;    // 144Xe
		case 1524: return 144.94472012;     // 145Xe
		case 1525: return 145.94851826;     // 146Xe
		case 1526: return 146.9542621;      // 147Xe
		case 1527: return 147.9581321;      // 148Xe
		case 1528: return 111.95030993;     // 112Cs
		case 1529: return 112.944429193;    // 113Cs
		case 1530: return 113.94129676;     // 114Cs
		case 1531: return 114.9359132;      // 115Cs
		case 1532: return 115.9333711;      // 116Cs
		case 1533: return 116.92861767;     // 117Cs
		case 1534: return 117.92656014;     // 118Cs
		case 1535: return 118.92237715;     // 119Cs
		case 1536: return 119.92067711;     // 120Cs
		case 1537: return 120.91722715;     // 121Cs
		case 1538: return 121.91610836;     // 122Cs
		case 1539: return 122.91299613;     // 123Cs
		case 1540: return 123.912257889;    // 124Cs
		case 1541: return 124.909728083;    // 125Cs
		case 1542: return 125.90944611;     // 126Cs
		case 1543: return 126.907417460;    // 127Cs
		case 1544: return 127.907748758;    // 128Cs
		case 1545: return 128.906065749;    // 129Cs
		case 1546: return 129.906709390;    // 130Cs
		case 1547: return 130.905464953;    // 131Cs
		case 1548: return 131.906433921;    // 132Cs
		case 1549: return 132.905451961080; // 133Cs
		case 1550: return 133.90671850317;  // 134Cs
		case 1551: return 134.905977011;    // 135Cs
		case 1552: return 135.907311420;    // 136Cs
		case 1553: return 136.9070892336;   // 137Cs
		case 1554: return 137.911017198;    // 138Cs
		case 1555: return 138.913363834;    // 139Cs
		case 1556: return 139.917283188;    // 140Cs
		case 1557: return 140.920045598;    // 141Cs
		case 1558: return 141.924296079;    // 142Cs
		case 1559: return 142.92734924;     // 143Cs
		case 1560: return 143.93207627;     // 144Cs
		case 1561: return 144.93552712;     // 145Cs
		case 1562: return 145.94034442;     // 146Cs
		case 1563: return 146.94415657;     // 147Cs
		case 1564: return 147.9492362;      // 148Cs
		case 1565: return 148.9530221;      // 149Cs
		case 1566: return 149.9583332;      // 150Cs
		case 1567: return 150.9625843;      // 151Cs
		case 1568: return 113.9506612;      // 114Ba
		case 1569: return 114.9473754;      // 115Ba
		case 1570: return 115.9412832;      // 116Ba
		case 1571: return 116.9381421;      // 117Ba
		case 1572: return 117.9330621;      // 118Ba
		case 1573: return 118.9306621;      // 119Ba
		case 1574: return 119.9260532;      // 120Ba
		case 1575: return 120.9240515;      // 121Ba
		case 1576: return 121.91990430;     // 122Ba
		case 1577: return 122.91878113;     // 123Ba
		case 1578: return 123.91509413;     // 124Ba
		case 1579: return 124.91447212;     // 125Ba
		case 1580: return 125.91125013;     // 126Ba
		case 1581: return 126.91109112;     // 127Ba
		case 1582: return 127.908342056;    // 128Ba
		case 1583: return 128.90868111;     // 129Ba
		case 1584: return 129.906320728;    // 130Ba
		case 1585: return 130.906941028;    // 131Ba
		case 1586: return 131.905061111;    // 132Ba
		case 1587: return 132.906007411;    // 133Ba
		case 1588: return 133.9045081830;   // 134Ba
		case 1589: return 134.9056883829;   // 135Ba
		case 1590: return 135.9045757329;   // 136Ba
		case 1591: return 136.9058271430;   // 137Ba
		case 1592: return 137.9052470031;   // 138Ba
		case 1593: return 138.9088411031;   // 139Ba
		case 1594: return 139.910605785;    // 140Ba
		case 1595: return 140.914403357;    // 141Ba
		case 1596: return 141.916432464;    // 142Ba
		case 1597: return 142.920625374;    // 143Ba
		case 1598: return 143.922954977;    // 144Ba
		case 1599: return 144.927518491;    // 145Ba
		case 1600: return 145.93028422;     // 146Ba
		case 1601: return 146.93530421;     // 147Ba
		case 1602: return 147.93817168;     // 148Ba
		case 1603: return 148.9430821;      // 149Ba
		case 1604: return 149.9460532;      // 150Ba
		case 1605: return 150.9512732;      // 151Ba
		case 1606: return 151.9548143;      // 152Ba
		case 1607: return 152.9603643;      // 153Ba
		case 1608: return 115.9563023;      // 116La
		case 1609: return 116.9499932;      // 117La
		case 1610: return 117.9467332;      // 118La
		case 1611: return 118.9409932;      // 119La
		case 1612: return 119.9380732;      // 120La
		case 1613: return 120.9331532;      // 121La
		case 1614: return 121.9307132;      // 122La
		case 1615: return 122.9263021;      // 123La
		case 1616: return 123.92457461;     // 124La
		case 1617: return 124.92081628;     // 125La
		case 1618: return 125.91951397;     // 126La
		case 1619: return 126.91637528;     // 127La
		case 1620: return 127.91559258;     // 128La
		case 1621: return 128.91269423;     // 129La
		case 1622: return 129.91236928;     // 130La
		case 1623: return 130.91007030;     // 131La
		case 1624: return 131.91011939;     // 132La
		case 1625: return 132.90821830;     // 133La
		case 1626: return 133.90851421;     // 134La
		case 1627: return 134.90698410;     // 135La
		case 1628: return 135.90763557;     // 136La
		case 1629: return 136.906450418;    // 137La
		case 1630: return 137.907114937;    // 138La
		case 1631: return 138.906356324;    // 139La
		case 1632: return 139.909480624;    // 140La
		case 1633: return 140.910966048;    // 141La
		case 1634: return 141.914090969;    // 142La
		case 1635: return 142.916079579;    // 143La
		case 1636: return 143.91964614;     // 144La
		case 1637: return 144.92180813;     // 145La
		case 1638: return 145.92587536;     // 146La
		case 1639: return 146.92841812;     // 147La
		case 1640: return 147.93267921;     // 148La
		case 1641: return 148.9353521;      // 149La
		case 1642: return 149.9394721;      // 150La
		case 1643: return 150.9423221;      // 151La
		case 1644: return 151.9468232;      // 152La
		case 1645: return 152.9503632;      // 153La
		case 1646: return 153.9551743;      // 154La
		case 1647: return 154.9590143;      // 155La
		case 1648: return 118.9527154;      // 119Ce
		case 1649: return 119.9465454;      // 120Ce
		case 1650: return 120.9433543;      // 121Ce
		case 1651: return 121.9378743;      // 122Ce
		case 1652: return 122.9352832;      // 123Ce
		case 1653: return 123.9303132;      // 124Ce
		case 1654: return 124.9284421;      // 125Ce
		case 1655: return 125.92397130;     // 126Ce
		case 1656: return 126.92272731;     // 127Ce
		case 1657: return 127.91891130;     // 128Ce
		case 1658: return 128.91810230;     // 129Ce
		case 1659: return 129.91473630;     // 130Ce
		case 1660: return 130.91442935;     // 131Ce
		case 1661: return 131.91146422;     // 132Ce
		case 1662: return 132.91152018;     // 133Ce
		case 1663: return 133.90892822;     // 134Ce
		case 1664: return 134.90916111;     // 135Ce
		case 1665: return 135.9071292141;   // 136Ce
		case 1666: return 136.9077623645;   // 137Ce
		case 1667: return 137.90599111;     // 138Ce
		case 1668: return 138.906655178;    // 139Ce
		case 1669: return 139.905443123;    // 140Ce
		case 1670: return 140.908280723;    // 141Ce
		case 1671: return 141.909250429;    // 142Ce
		case 1672: return 142.912392129;    // 143Ce
		case 1673: return 143.913652934;    // 144Ce
		case 1674: return 144.91726536;     // 145Ce
		case 1675: return 145.91880218;     // 146Ce
		case 1676: return 146.922689992;    // 147Ce
		case 1677: return 147.92442412;     // 148Ce
		case 1678: return 148.92842711;     // 149Ce
		case 1679: return 149.93038413;     // 150Ce
		case 1680: return 150.93427219;     // 151Ce
		case 1681: return 151.9366021;      // 152Ce
		case 1682: return 152.9409321;      // 153Ce
		case 1683: return 153.9438032;      // 154Ce
		case 1684: return 154.9485543;      // 155Ce
		case 1685: return 155.9518343;      // 156Ce
		case 1686: return 156.9570554;      // 157Ce
		case 1687: return 120.9553254;      // 121Pr
		case 1688: return 121.9517554;      // 122Pr
		case 1689: return 122.9459643;      // 123Pr
		case 1690: return 123.9429443;      // 124Pr
		case 1691: return 124.9377032;      // 125Pr
		case 1692: return 125.9352421;      // 126Pr
		case 1693: return 126.9307121;      // 127Pr
		case 1694: return 127.92879132;     // 128Pr
		case 1695: return 128.92509532;     // 129Pr
		case 1696: return 129.92359069;     // 130Pr
		case 1697: return 130.92023550;     // 131Pr
		case 1698: return 131.91925561;     // 132Pr
		case 1699: return 132.91633113;     // 133Pr
		case 1700: return 133.91569722;     // 134Pr
		case 1701: return 134.91311213;     // 135Pr
		case 1702: return 135.91267712;     // 136Pr
		case 1703: return 136.910679287;    // 137Pr
		case 1704: return 137.91075415;     // 138Pr
		case 1705: return 138.908940885;    // 139Pr
		case 1706: return 139.909080369;    // 140Pr
		case 1707: return 140.907657623;    // 141Pr
		case 1708: return 141.910049623;    // 142Pr
		case 1709: return 142.910822824;    // 143Pr
		case 1710: return 143.913310932;    // 144Pr
		case 1711: return 144.914518278;    // 145Pr
		case 1712: return 145.91768037;     // 146Pr
		case 1713: return 146.91900817;     // 147Pr
		case 1714: return 147.92213016;     // 148Pr
		case 1715: return 148.92373611;     // 149Pr
		case 1716: return 149.926676597;    // 150Pr
		case 1717: return 150.92830913;     // 151Pr
		case 1718: return 151.93155320;     // 152Pr
		case 1719: return 152.93390413;     // 153Pr
		case 1720: return 153.9375316;      // 154Pr
		case 1721: return 154.94050918;     // 155Pr
		case 1722: return 155.9446432;      // 156Pr
		case 1723: return 156.9478943;      // 157Pr
		case 1724: return 157.9524143;      // 158Pr
		case 1725: return 158.9558954;      // 159Pr
		case 1726: return 123.9522054;      // 124Nd
		case 1727: return 124.9489043;      // 125Nd
		case 1728: return 125.9431132;      // 126Nd
		case 1729: return 126.9403832;      // 127Nd
		case 1730: return 127.9352521;      // 128Nd
		case 1731: return 128.9331022;      // 129Nd
		case 1732: return 129.92850630;     // 130Nd
		case 1733: return 130.92724830;     // 131Nd
		case 1734: return 131.92332126;     // 132Nd
		case 1735: return 132.92234850;     // 133Nd
		case 1736: return 133.91879013;     // 134Nd
		case 1737: return 134.91818121;     // 135Nd
		case 1738: return 135.91497613;     // 136Nd
		case 1739: return 136.91456213;     // 137Nd
		case 1740: return 137.91195012;     // 138Nd
		case 1741: return 138.91195430;     // 139Nd
		case 1742: return 139.90955028;     // 140Nd
		case 1743: return 140.909614738;    // 141Nd
		case 1744: return 141.907729020;    // 142Nd
		case 1745: return 142.909820020;    // 143Nd
		case 1746: return 143.910093020;    // 144Nd
		case 1747: return 144.912579320;    // 145Nd
		case 1748: return 145.913122620;    // 146Nd
		case 1749: return 146.916106120;    // 147Nd
		case 1750: return 147.916899326;    // 148Nd
		case 1751: return 148.920154826;    // 149Nd
		case 1752: return 149.920902218;    // 150Nd
		case 1753: return 150.923840318;    // 151Nd
		case 1754: return 151.92469226;     // 152Nd
		case 1755: return 152.927718029;    // 153Nd
		case 1756: return 153.9294812;      // 154Nd
		case 1757: return 154.933135798;    // 155Nd
		case 1758: return 155.9350821;      // 156Nd
		case 1759: return 156.93938627;     // 157Nd
		case 1760: return 157.9419732;      // 158Nd
		case 1761: return 158.9465343;      // 159Nd
		case 1762: return 159.9494043;      // 160Nd
		case 1763: return 160.9542854;      // 161Nd
		case 1764: return 125.9579254;      // 126Pm
		case 1765: return 126.9519243;      // 127Pm
		case 1766: return 127.9487032;      // 128Pm
		case 1767: return 128.9432332;      // 129Pm
		case 1768: return 129.9405321;      // 130Pm
		case 1769: return 130.9356721;      // 131Pm
		case 1770: return 131.9338416;      // 132Pm
		case 1771: return 132.92978254;     // 133Pm
		case 1772: return 133.92835362;     // 134Pm
		case 1773: return 134.92482370;     // 135Pm
		case 1774: return 135.92358577;     // 136Pm
		case 1775: return 136.92048014;     // 137Pm
		case 1776: return 137.91954830;     // 138Pm
		case 1777: return 138.91680015;     // 139Pm
		case 1778: return 139.91604038;     // 140Pm
		case 1779: return 140.91355515;     // 141Pm
		case 1780: return 141.91289025;     // 142Pm
		case 1781: return 142.910938334;    // 143Pm
		case 1782: return 143.912596434;    // 144Pm
		case 1783: return 144.912755933;    // 145Pm
		case 1784: return 145.914702448;    // 146Pm
		case 1785: return 146.915145019;    // 147Pm
		case 1786: return 147.917481963;    // 148Pm
		case 1787: return 148.918342327;    // 149Pm
		case 1788: return 149.92099122;     // 150Pm
		case 1789: return 150.921217551;    // 151Pm
		case 1790: return 151.92350628;     // 152Pm
		case 1791: return 152.924156797;    // 153Pm
		case 1792: return 153.92647248;     // 154Pm
		case 1793: return 154.928137051;    // 155Pm
		case 1794: return 155.931117539;    // 156Pm
		case 1795: return 156.933121475;    // 157Pm
		case 1796: return 157.93656514;     // 158Pm
		case 1797: return 158.93928711;     // 159Pm
		case 1798: return 159.9431032;      // 160Pm
		case 1799: return 160.9460732;      // 161Pm
		case 1800: return 161.9502243;      // 162Pm
		case 1801: return 162.9535754;      // 163Pm
		case 1802: return 127.9584254;      // 128Sm
		case 1803: return 128.9547654;      // 129Sm
		case 1804: return 129.9490043;      // 130Sm
		case 1805: return 130.9461843;      // 131Sm
		case 1806: return 131.9408732;      // 132Sm
		case 1807: return 132.9385632;      // 133Sm
		case 1808: return 133.9341121;      // 134Sm
		case 1809: return 134.9325217;      // 135Sm
		case 1810: return 135.92827613;     // 136Sm
		case 1811: return 136.92697146;     // 137Sm
		case 1812: return 137.92324413;     // 138Sm
		case 1813: return 138.92229712;     // 139Sm
		case 1814: return 139.91899513;     // 140Sm
		case 1815: return 140.918481692;    // 141Sm
		case 1816: return 141.915204436;    // 142Sm
		case 1817: return 142.914635333;    // 143Sm
		case 1818: return 143.912006521;    // 144Sm
		case 1819: return 144.913417321;    // 145Sm
		case 1820: return 145.913047035;    // 146Sm
		case 1821: return 146.914904419;    // 147Sm
		case 1822: return 147.914829219;    // 148Sm
		case 1823: return 148.917192118;    // 149Sm
		case 1824: return 149.917282918;    // 150Sm
		case 1825: return 150.919939818;    // 151Sm
		case 1826: return 151.919739718;    // 152Sm
		case 1827: return 152.922104718;    // 153Sm
		case 1828: return 153.922216920;    // 154Sm
		case 1829: return 154.924647720;    // 155Sm
		case 1830: return 155.92553610;     // 156Sm
		case 1831: return 156.928418748;    // 157Sm
		case 1832: return 157.929951053;    // 158Sm
		case 1833: return 158.933217264;    // 159Sm
		case 1834: return 159.935335364;    // 160Sm
		case 1835: return 160.939160273;    // 161Sm
		case 1836: return 161.9414621;      // 162Sm
		case 1837: return 162.9455532;      // 163Sm
		case 1838: return 163.9483632;      // 164Sm
		case 1839: return 164.9529743;      // 165Sm
		case 1840: return 129.9636954;      // 130Eu
		case 1841: return 130.9578443;      // 131Eu
		case 1842: return 131.9546743;      // 132Eu
		case 1843: return 132.9492932;      // 133Eu
		case 1844: return 133.9464032;      // 134Eu
		case 1845: return 134.9418721;      // 135Eu
		case 1846: return 135.9396221;      // 136Eu
		case 1847: return 136.9354621;      // 137Eu
		case 1848: return 137.93370930;     // 138Eu
		case 1849: return 138.92979214;     // 139Eu
		case 1850: return 139.92808855;     // 140Eu
		case 1851: return 140.92493214;     // 141Eu
		case 1852: return 141.92344232;     // 142Eu
		case 1853: return 142.92029912;     // 143Eu
		case 1854: return 143.91882012;     // 144Eu
		case 1855: return 144.916272636;    // 145Eu
		case 1856: return 145.917211065;    // 146Eu
		case 1857: return 146.916752731;    // 147Eu
		case 1858: return 147.91808911;     // 148Eu
		case 1859: return 148.917937844;    // 149Eu
		case 1860: return 149.919707768;    // 150Eu
		case 1861: return 150.919857818;    // 151Eu
		case 1862: return 151.921752218;    // 152Eu
		case 1863: return 152.921238018;    // 153Eu
		case 1864: return 153.922987019;    // 154Eu
		case 1865: return 154.922901119;    // 155Eu
		case 1866: return 155.924760559;    // 156Eu
		case 1867: return 156.925433446;    // 157Eu
		case 1868: return 157.92779911;     // 158Eu
		case 1869: return 158.929100147;    // 159Eu
		case 1870: return 159.93185110;     // 160Eu
		case 1871: return 160.93366411;     // 161Eu
		case 1872: return 161.93698965;     // 162Eu
		case 1873: return 162.93919676;     // 163Eu
		case 1874: return 163.9427422;      // 164Eu
		case 1875: return 164.9455935;      // 165Eu
		case 1876: return 165.9496232;      // 166Eu
		case 1877: return 166.9528943;      // 167Eu
		case 1878: return 132.9613354;      // 133Gd
		case 1879: return 133.9556643;      // 134Gd
		case 1880: return 134.9524543;      // 135Gd
		case 1881: return 135.9473032;      // 136Gd
		case 1882: return 136.9450232;      // 137Gd
		case 1883: return 137.9402521;      // 138Gd
		case 1884: return 138.9381321;      // 139Gd
		case 1885: return 139.93367430;     // 140Gd
		case 1886: return 140.93212621;     // 141Gd
		case 1887: return 141.92811630;     // 142Gd
		case 1888: return 142.9267522;      // 143Gd
		case 1889: return 143.92296330;     // 144Gd
		case 1890: return 144.92171321;     // 145Gd
		case 1891: return 145.918318846;    // 146Gd
		case 1892: return 146.919101425;    // 147Gd
		case 1893: return 147.918121521;    // 148Gd
		case 1894: return 148.919348138;    // 149Gd
		case 1895: return 149.918664466;    // 150Gd
		case 1896: return 150.920356035;    // 151Gd
		case 1897: return 151.919799518;    // 152Gd
		case 1898: return 152.921758018;    // 153Gd
		case 1899: return 153.920874117;    // 154Gd
		case 1900: return 154.922630517;    // 155Gd
		case 1901: return 155.922131217;    // 156Gd
		case 1902: return 156.923968617;    // 157Gd
		case 1903: return 157.924112317;    // 158Gd
		case 1904: return 158.926397017;    // 159Gd
		case 1905: return 159.927062418;    // 160Gd
		case 1906: return 160.929677521;    // 161Gd
		case 1907: return 161.930993045;    // 162Gd
		case 1908: return 162.934176990;    // 163Gd
		case 1909: return 163.9358321;      // 164Gd
		case 1910: return 164.9393632;      // 165Gd
		case 1911: return 165.9414664;      // 166Gd
		case 1912: return 166.9454543;      // 167Gd
		case 1913: return 167.9480843;      // 168Gd
		case 1914: return 168.9526054;      // 169Gd
		case 1915: return 134.9647643;      // 135Tb
		case 1916: return 135.9612954;      // 136Tb
		case 1917: return 136.9560254;      // 137Tb
		case 1918: return 137.9531232;      // 138Tb
		case 1919: return 138.9483332;      // 139Tb
		case 1920: return 139.9458186;      // 140Tb
		case 1921: return 140.9414511;      // 141Tb
		case 1922: return 141.9392875;      // 142Tb
		case 1923: return 142.93513755;     // 143Tb
		case 1924: return 143.93304530;     // 144Tb
		case 1925: return 144.9288210;      // 145Tb
		case 1926: return 145.92725348;     // 146Tb
		case 1927: return 146.924054887;    // 147Tb
		case 1928: return 147.92428214;     // 148Tb
		case 1929: return 148.923253541;    // 149Tb
		case 1930: return 149.923664980;    // 150Tb
		case 1931: return 150.923109646;    // 151Tb
		case 1932: return 151.92408343;     // 152Tb
		case 1933: return 152.923442444;    // 153Tb
		case 1934: return 153.92468549;     // 154Tb
		case 1935: return 154.92351111;     // 155Tb
		case 1936: return 155.924755243;    // 156Tb
		case 1937: return 156.924033018;    // 157Tb
		case 1938: return 157.925420920;    // 158Tb
		case 1939: return 158.925354719;    // 159Tb
		case 1940: return 159.927175619;    // 160Tb
		case 1941: return 160.927577820;    // 161Tb
		case 1942: return 161.92949539;     // 162Tb
		case 1943: return 162.930654747;    // 163Tb
		case 1944: return 163.9333611;      // 164Tb
		case 1945: return 164.9349821;      // 165Tb
		case 1946: return 165.93786075;     // 166Tb
		case 1947: return 166.9399621;      // 167Tb
		case 1948: return 167.9434032;      // 168Tb
		case 1949: return 168.9459732;      // 169Tb
		case 1950: return 169.9498443;      // 170Tb
		case 1951: return 170.9527354;      // 171Tb
		case 1952: return 137.9625043;      // 138Dy
		case 1953: return 138.9595954;      // 139Dy
		case 1954: return 139.9540254;      // 140Dy
		case 1955: return 140.9512832;      // 141Dy
		case 1956: return 141.9461978;      // 142Dy
		case 1957: return 142.94399414;     // 143Dy
		case 1958: return 143.939269577;    // 144Dy
		case 1959: return 144.937474070;    // 145Dy
		case 1960: return 145.932844572;    // 146Dy
		case 1961: return 146.931082795;    // 147Dy
		case 1962: return 147.92715710;     // 148Dy
		case 1963: return 148.92732210;     // 149Dy
		case 1964: return 149.925593348;    // 150Dy
		case 1965: return 150.926191638;    // 151Dy
		case 1966: return 151.924725351;    // 152Dy
		case 1967: return 152.925772445;    // 153Dy
		case 1968: return 153.924429380;    // 154Dy
		case 1969: return 154.92575910;     // 155Dy
		case 1970: return 155.924284717;    // 156Dy
		case 1971: return 156.925470757;    // 157Dy
		case 1972: return 157.924415931;    // 158Dy
		case 1973: return 158.925747022;    // 159Dy
		case 1974: return 159.925204620;    // 160Dy
		case 1975: return 160.926940520;    // 161Dy
		case 1976: return 161.926805620;    // 162Dy
		case 1977: return 162.928738320;    // 163Dy
		case 1978: return 163.929181920;    // 164Dy
		case 1979: return 164.931710520;    // 165Dy
		case 1980: return 165.932813921;    // 166Dy
		case 1981: return 166.93566165;     // 167Dy
		case 1982: return 167.9371315;      // 168Dy
		case 1983: return 168.9403132;      // 169Dy
		case 1984: return 169.9423921;      // 170Dy
		case 1985: return 170.9461232;      // 171Dy
		case 1986: return 171.9484632;      // 172Dy
		case 1987: return 172.9528343;      // 173Dy
		case 1988: return 139.9685954;      // 140Ho
		case 1989: return 140.9631154;      // 141Ho
		case 1990: return 141.9600154;      // 142Ho
		case 1991: return 142.9548643;      // 143Ho
		case 1992: return 143.952109791;    // 144Ho
		case 1993: return 144.947267480;    // 145Ho
		case 1994: return 145.944993571;    // 146Ho
		case 1995: return 146.940142354;    // 147Ho
		case 1996: return 147.93774490;     // 148Ho
		case 1997: return 148.93380316;     // 149Ho
		case 1998: return 149.93349815;     // 150Ho
		case 1999: return 150.931698389;    // 151Ho
		case 2000: return 151.93172414;     // 152Ho
		case 2001: return 152.930206456;    // 153Ho
		case 2002: return 153.930606889;    // 154Ho
		case 2003: return 154.92910419;     // 155Ho
		case 2004: return 155.92970664;     // 156Ho
		case 2005: return 156.92825425;     // 157Ho
		case 2006: return 157.92894629;     // 158Ho
		case 2007: return 158.927719736;    // 159Ho
		case 2008: return 159.92873716;     // 160Ho
		case 2009: return 160.927861530;    // 161Ho
		case 2010: return 161.929102339;    // 162Ho
		case 2011: return 162.928741020;    // 163Ho
		case 2012: return 163.930240325;    // 164Ho
		case 2013: return 164.930328821;    // 165Ho
		case 2014: return 165.932290921;    // 166Ho
		case 2015: return 166.933138559;    // 167Ho
		case 2016: return 167.93552232;     // 168Ho
		case 2017: return 168.93687822;     // 169Ho
		case 2018: return 169.93962554;     // 170Ho
		case 2019: return 170.9414764;      // 171Ho
		case 2020: return 171.9447321;      // 172Ho
		case 2021: return 172.9470232;      // 173Ho
		case 2022: return 173.9509532;      // 174Ho
		case 2023: return 174.9536243;      // 175Ho
		case 2024: return 141.9701054;      // 142Er
		case 2025: return 142.9666243;      // 143Er
		case 2026: return 143.9607021;      // 144Er
		case 2027: return 144.9580521;      // 145Er
		case 2028: return 145.952418472;    // 146Er
		case 2029: return 146.94996441;     // 147Er
		case 2030: return 147.94473511;     // 148Er
		case 2031: return 148.94230630;     // 149Er
		case 2032: return 149.93791618;     // 150Er
		case 2033: return 150.93744918;     // 151Er
		case 2034: return 151.93505710;     // 152Er
		case 2035: return 152.93508010;     // 153Er
		case 2036: return 153.932790855;    // 154Er
		case 2037: return 154.933215967;    // 155Er
		case 2038: return 155.93106726;     // 156Er
		case 2039: return 156.93194927;     // 157Er
		case 2040: return 157.92989327;     // 158Er
		case 2041: return 158.930691842;    // 159Er
		case 2042: return 159.92907726;     // 160Er
		case 2043: return 160.930004696;    // 161Er
		case 2044: return 161.928788420;    // 162Er
		case 2045: return 162.930040853;    // 163Er
		case 2046: return 163.929208820;    // 164Er
		case 2047: return 164.930734521;    // 165Er
		case 2048: return 165.930299522;    // 166Er
		case 2049: return 166.932054622;    // 167Er
		case 2050: return 167.932376722;    // 168Er
		case 2051: return 168.934596822;    // 169Er
		case 2052: return 169.935470226;    // 170Er
		case 2053: return 170.938035726;    // 171Er
		case 2054: return 171.939361947;    // 172Er
		case 2055: return 172.9424021;      // 173Er
		case 2056: return 173.9442332;      // 174Er
		case 2057: return 174.9477743;      // 175Er
		case 2058: return 175.9499443;      // 176Er
		case 2059: return 176.9539954;      // 177Er
		case 2060: return 143.9762843;      // 144Tm
		case 2061: return 144.9703921;      // 145Tm
		case 2062: return 145.9668421;      // 146Tm
		case 2063: return 146.961379973;    // 147Tm
		case 2064: return 147.95838411;     // 148Tm
		case 2065: return 148.9528932;      // 149Tm
		case 2066: return 149.9500921;      // 150Tm
		case 2067: return 150.94548821;     // 151Tm
		case 2068: return 151.94442279;     // 152Tm
		case 2069: return 152.94204016;     // 153Tm
		case 2070: return 153.94157015;     // 154Tm
		case 2071: return 154.93921011;     // 155Tm
		case 2072: return 155.93899216;     // 156Tm
		case 2073: return 156.93694428;     // 157Tm
		case 2074: return 157.93698027;     // 158Tm
		case 2075: return 158.93497530;     // 159Tm
		case 2076: return 159.93526337;     // 160Tm
		case 2077: return 160.93354930;     // 161Tm
		case 2078: return 161.93400228;     // 162Tm
		case 2079: return 162.932659262;    // 163Tm
		case 2080: return 163.93354426;     // 164Tm
		case 2081: return 164.932443126;    // 165Tm
		case 2082: return 165.93356113;     // 166Tm
		case 2083: return 166.932856225;    // 167Tm
		case 2084: return 167.934177427;    // 168Tm
		case 2085: return 168.934217922;    // 169Tm
		case 2086: return 169.935806022;    // 170Tm
		case 2087: return 170.936433924;    // 171Tm
		case 2088: return 171.938405562;    // 172Tm
		case 2089: return 172.939608453;    // 173Tm
		case 2090: return 173.94217348;     // 174Tm
		case 2091: return 174.94384154;     // 175Tm
		case 2092: return 175.9470011;      // 176Tm
		case 2093: return 176.9490432;      // 177Tm
		case 2094: return 177.9526443;      // 178Tm
		case 2095: return 178.9553454;      // 179Tm
		case 2096: return 147.9675864;      // 148Yb
		case 2097: return 148.9643654;      // 149Yb
		case 2098: return 149.9585243;      // 150Yb
		case 2099: return 150.9554032;      // 151Yb
		case 2100: return 151.9502717;      // 152Yb
		case 2101: return 152.9493221;      // 153Yb
		case 2102: return 153.94639619;     // 154Yb
		case 2103: return 154.94578318;     // 155Yb
		case 2104: return 155.94282511;     // 156Yb
		case 2105: return 156.94264512;     // 157Yb
		case 2106: return 157.939870586;    // 158Yb
		case 2107: return 158.94005519;     // 159Yb
		case 2108: return 159.93755717;     // 160Yb
		case 2109: return 160.93790717;     // 161Yb
		case 2110: return 161.93577417;     // 162Yb
		case 2111: return 162.93634017;     // 163Yb
		case 2112: return 163.93449517;     // 164Yb
		case 2113: return 164.93527028;     // 165Yb
		case 2114: return 165.933874778;    // 166Yb
		case 2115: return 166.934953047;    // 167Yb
		case 2116: return 167.933889622;    // 168Yb
		case 2117: return 168.935182522;    // 169Yb
		case 2118: return 169.934766422;    // 170Yb
		case 2119: return 170.936330222;    // 171Yb
		case 2120: return 171.936385922;    // 172Yb
		case 2121: return 172.938215122;    // 173Yb
		case 2122: return 173.938866422;    // 174Yb
		case 2123: return 174.941280822;    // 175Yb
		case 2124: return 175.942576424;    // 176Yb
		case 2125: return 176.945265624;    // 177Yb
		case 2126: return 177.94665111;     // 178Yb
		case 2127: return 178.9500421;      // 179Yb
		case 2128: return 179.9521232;      // 180Yb
		case 2129: return 180.9558932;      // 181Yb
		case 2130: return 149.9735554;      // 150Lu
		case 2131: return 150.9676843;      // 151Lu
		case 2132: return 151.9641221;      // 152Lu
		case 2133: return 152.9587517;      // 153Lu
		case 2134: return 153.9573622;      // 154Lu
		case 2135: return 154.95432121;     // 155Lu
		case 2136: return 155.95303379;     // 156Lu
		case 2137: return 156.95012716;     // 157Lu
		case 2138: return 157.94931616;     // 158Lu
		case 2139: return 158.94663640;     // 159Lu
		case 2140: return 159.94603361;     // 160Lu
		case 2141: return 160.94357230;     // 161Lu
		case 2142: return 161.94328381;     // 162Lu
		case 2143: return 162.94117930;     // 163Lu
		case 2144: return 163.94133930;     // 164Lu
		case 2145: return 164.93940728;     // 165Lu
		case 2146: return 165.93985932;     // 166Lu
		case 2147: return 166.93827034;     // 167Lu
		case 2148: return 167.93873642;     // 168Lu
		case 2149: return 168.937644139;    // 169Lu
		case 2150: return 169.93847818;     // 170Lu
		case 2151: return 170.937917027;    // 171Lu
		case 2152: return 171.939089130;    // 172Lu
		case 2153: return 172.938934023;    // 173Lu
		case 2154: return 173.940340923;    // 174Lu
		case 2155: return 174.940775220;    // 175Lu
		case 2156: return 175.942689720;    // 176Lu
		case 2157: return 176.943761520;    // 177Lu
		case 2158: return 177.945958029;    // 178Lu
		case 2159: return 178.947330957;    // 179Lu
		case 2160: return 179.94988876;     // 180Lu
		case 2161: return 180.9519117;      // 181Lu
		case 2162: return 181.9550421;      // 182Lu
		case 2163: return 182.95736398;     // 183Lu
		case 2164: return 183.9609132;      // 184Lu
		case 2165: return 184.9636232;      // 185Lu
		case 2166: return 152.9706954;      // 153Hf
		case 2167: return 153.9648654;      // 154Hf
		case 2168: return 154.9631132;      // 155Hf
		case 2169: return 155.9593517;      // 156Hf
		case 2170: return 156.9582421;      // 157Hf
		case 2171: return 157.95480119;     // 158Hf
		case 2172: return 158.95399618;     // 159Hf
		case 2173: return 159.95069111;     // 160Hf
		case 2174: return 160.95027824;     // 161Hf
		case 2175: return 161.947214897;    // 162Hf
		case 2176: return 162.94711327;     // 163Hf
		case 2177: return 163.94437117;     // 164Hf
		case 2178: return 164.94456730;     // 165Hf
		case 2179: return 165.94218030;     // 166Hf
		case 2180: return 166.94260030;     // 167Hf
		case 2181: return 167.94056830;     // 168Hf
		case 2182: return 168.94125930;     // 169Hf
		case 2183: return 169.93960930;     // 170Hf
		case 2184: return 170.94049231;     // 171Hf
		case 2185: return 171.93945026;     // 172Hf
		case 2186: return 172.94051330;     // 173Hf
		case 2187: return 173.940046128;    // 174Hf
		case 2188: return 174.941509229;    // 175Hf
		case 2189: return 175.941407622;    // 176Hf
		case 2190: return 176.943227720;    // 177Hf
		case 2191: return 177.943705820;    // 178Hf
		case 2192: return 178.945823220;    // 179Hf
		case 2193: return 179.946557020;    // 180Hf
		case 2194: return 180.949108320;    // 181Hf
		case 2195: return 181.950561268;    // 182Hf
		case 2196: return 182.95353032;     // 183Hf
		case 2197: return 183.95544643;     // 184Hf
		case 2198: return 184.95886298;     // 185Hf
		case 2199: return 185.96089759;     // 186Hf
		case 2200: return 186.9647732;      // 187Hf
		case 2201: return 187.9668532;      // 188Hf
		case 2202: return 188.9708432;      // 189Hf
		case 2203: return 154.9742454;      // 155Ta
		case 2204: return 155.9720332;      // 156Ta
		case 2205: return 156.9681817;      // 157Ta
		case 2206: return 157.9665422;      // 158Ta
		case 2207: return 158.96302321;     // 159Ta
		case 2208: return 159.96148879;     // 160Ta
		case 2209: return 160.95845227;     // 161Ta
		case 2210: return 161.95729456;     // 162Ta
		case 2211: return 162.95433741;     // 163Ta
		case 2212: return 163.95353430;     // 164Ta
		case 2213: return 164.95078115;     // 165Ta
		case 2214: return 165.95051230;     // 166Ta
		case 2215: return 166.94809330;     // 167Ta
		case 2216: return 167.94804730;     // 168Ta
		case 2217: return 168.94601130;     // 169Ta
		case 2218: return 169.94617530;     // 170Ta
		case 2219: return 170.94447630;     // 171Ta
		case 2220: return 171.94489530;     // 172Ta
		case 2221: return 172.94375030;     // 173Ta
		case 2222: return 173.94445430;     // 174Ta
		case 2223: return 174.94373730;     // 175Ta
		case 2224: return 175.94485733;     // 176Ta
		case 2225: return 176.944479538;    // 177Ta
		case 2226: return 177.94567856;     // 178Ta
		case 2227: return 178.945936621;    // 179Ta
		case 2228: return 179.947464824;    // 180Ta
		case 2229: return 180.947995820;    // 181Ta
		case 2230: return 181.950151920;    // 182Ta
		case 2231: return 182.951372620;    // 183Ta
		case 2232: return 183.95400828;     // 184Ta
		case 2233: return 184.95555915;     // 185Ta
		case 2234: return 185.95855164;     // 186Ta
		case 2235: return 186.96038671;     // 187Ta
		case 2236: return 187.96391671;     // 188Ta
		case 2237: return 188.9658332;      // 189Ta
		case 2238: return 189.9693921;      // 190Ta
		case 2239: return 190.9715632;      // 191Ta
		case 2240: return 191.9751443;      // 192Ta
		case 2241: return 156.9788443;      // 157W
		case 2242: return 157.9745654;      // 158W
		case 2243: return 158.9726432;      // 159W
		case 2244: return 159.9684617;      // 160W
		case 2245: return 160.9672021;      // 161W
		case 2246: return 161.96349919;     // 162W
		case 2247: return 162.96252457;     // 163W
		case 2248: return 163.95896111;     // 164W
		case 2249: return 164.95828127;     // 165W
		case 2250: return 165.95503110;     // 166W
		case 2251: return 166.95480520;     // 167W
		case 2252: return 167.95180614;     // 168W
		case 2253: return 168.95177917;     // 169W
		case 2254: return 169.94923214;     // 170W
		case 2255: return 170.94945130;     // 171W
		case 2256: return 171.94729230;     // 172W
		case 2257: return 172.94768930;     // 173W
		case 2258: return 173.94607930;     // 174W
		case 2259: return 174.94671730;     // 175W
		case 2260: return 175.94563430;     // 176W
		case 2261: return 176.94664330;     // 177W
		case 2262: return 177.94588316;     // 178W
		case 2263: return 178.94707716;     // 179W
		case 2264: return 179.946710820;    // 180W
		case 2265: return 180.948197851;    // 181W
		case 2266: return 181.9482039491;   // 182W
		case 2267: return 182.9502227590;   // 183W
		case 2268: return 183.9509309294;   // 184W
		case 2269: return 184.9534189799;   // 185W
		case 2270: return 185.954362817;    // 186W
		case 2271: return 186.957158817;    // 187W
		case 2272: return 187.958486236;    // 188W
		case 2273: return 188.96176344;     // 189W
		case 2274: return 189.96309142;     // 190W
		case 2275: return 190.96653148;     // 191W
		case 2276: return 191.9681721;      // 192W
		case 2277: return 192.9717821;      // 193W
		case 2278: return 193.9736732;      // 194W
		case 2279: return 158.9841854;      // 159Re
		case 2280: return 159.9818232;      // 160Re
		case 2281: return 160.9775717;      // 161Re
		case 2282: return 161.9758422;      // 162Re
		case 2283: return 162.97208020;     // 163Re
		case 2284: return 163.97045379;     // 164Re
		case 2285: return 164.96710327;     // 165Re
		case 2286: return 165.96576178;     // 166Re
		case 2287: return 166.96259544;     // 167Re
		case 2288: return 167.96157333;     // 168Re
		case 2289: return 168.95876612;     // 169Re
		case 2290: return 169.95822028;     // 170Re
		case 2291: return 170.95571630;     // 171Re
		case 2292: return 171.95542042;     // 172Re
		case 2293: return 172.95324330;     // 173Re
		case 2294: return 173.95311530;     // 174Re
		case 2295: return 174.95138130;     // 175Re
		case 2296: return 175.95162330;     // 176Re
		case 2297: return 176.95032830;     // 177Re
		case 2298: return 177.95098930;     // 178Re
		case 2299: return 178.94998926;     // 179Re
		case 2300: return 179.95079223;     // 180Re
		case 2301: return 180.95005814;     // 181Re
		case 2302: return 181.9512111;      // 182Re
		case 2303: return 182.950819686;    // 183Re
		case 2304: return 183.952522847;    // 184Re
		case 2305: return 184.952954513;    // 185Re
		case 2306: return 185.954985613;    // 186Re
		case 2307: return 186.955750116;    // 187Re
		case 2308: return 187.958111516;    // 188Re
		case 2309: return 188.959226089;    // 189Re
		case 2310: return 189.96174476;     // 190Re
		case 2311: return 190.96312211;     // 191Re
		case 2312: return 191.96608882;     // 192Re
		case 2313: return 192.96754141;     // 193Re
		case 2314: return 193.9707621;      // 194Re
		case 2315: return 194.9725432;      // 195Re
		case 2316: return 195.9758032;      // 196Re
		case 2317: return 196.9779932;      // 197Re
		case 2318: return 197.9816043;      // 198Re
		case 2319: return 160.9890343;      // 161Os
		case 2320: return 161.9844354;      // 162Os
		case 2321: return 162.9824132;      // 163Os
		case 2322: return 163.9780217;      // 164Os
		case 2323: return 164.9766022;      // 165Os
		case 2324: return 165.97269220;     // 166Os
		case 2325: return 166.97154978;     // 167Os
		case 2326: return 167.96780812;     // 168Os
		case 2327: return 168.96701827;     // 169Os
		case 2328: return 169.96357811;     // 170Os
		case 2329: return 170.96317419;     // 171Os
		case 2330: return 171.96001714;     // 172Os
		case 2331: return 172.95980816;     // 173Os
		case 2332: return 173.95706411;     // 174Os
		case 2333: return 174.95694513;     // 175Os
		case 2334: return 175.95480630;     // 176Os
		case 2335: return 176.95496617;     // 177Os
		case 2336: return 177.95325415;     // 178Os
		case 2337: return 178.95381718;     // 179Os
		case 2338: return 179.95237517;     // 180Os
		case 2339: return 180.95324727;     // 181Os
		case 2340: return 181.95211023;     // 182Os
		case 2341: return 182.95312553;     // 183Os
		case 2342: return 183.952488514;    // 184Os
		case 2343: return 184.954041714;    // 185Os
		case 2344: return 185.953835016;    // 186Os
		case 2345: return 186.955747416;    // 187Os
		case 2346: return 187.955835216;    // 188Os
		case 2347: return 188.958144217;    // 189Os
		case 2348: return 189.958443717;    // 190Os
		case 2349: return 190.960926417;    // 191Os
		case 2350: return 191.961477029;    // 192Os
		case 2351: return 192.964147929;    // 193Os
		case 2352: return 193.965177230;    // 194Os
		case 2353: return 194.96831865;     // 195Os
		case 2354: return 195.96964143;     // 196Os
		case 2355: return 196.9728321;      // 197Os
		case 2356: return 197.9744121;      // 198Os
		case 2357: return 198.9780121;      // 199Os
		case 2358: return 199.9798432;      // 200Os
		case 2359: return 200.9836432;      // 201Os
		case 2360: return 201.9859543;      // 202Os
		case 2361: return 163.9919134;      // 164Ir
		case 2362: return 164.9875018;      // 165Ir
		case 2363: return 165.9856622;      // 166Ir
		case 2364: return 166.98166620;     // 167Ir
		case 2365: return 167.97990780;     // 168Ir
		case 2366: return 168.97629827;     // 169Ir
		case 2367: return 169.97492295;     // 170Ir
		case 2368: return 170.97164042;     // 171Ir
		case 2369: return 171.97060735;     // 172Ir
		case 2370: return 172.96750612;     // 173Ir
		case 2371: return 173.96686130;     // 174Ir
		case 2372: return 174.96415013;     // 175Ir
		case 2373: return 175.96365022;     // 176Ir
		case 2374: return 176.96130121;     // 177Ir
		case 2375: return 177.96108221;     // 178Ir
		case 2376: return 178.95912010;     // 179Ir
		case 2377: return 179.95922923;     // 180Ir
		case 2378: return 180.95762528;     // 181Ir
		case 2379: return 181.95807623;     // 182Ir
		case 2380: return 182.95684026;     // 183Ir
		case 2381: return 183.95747630;     // 184Ir
		case 2382: return 184.95669830;     // 185Ir
		case 2383: return 185.95794418;     // 186Ir
		case 2384: return 186.95754230;     // 187Ir
		case 2385: return 187.95882810;     // 188Ir
		case 2386: return 188.95871514;     // 189Ir
		case 2387: return 189.960541221;    // 190Ir
		case 2388: return 190.960589321;    // 191Ir
		case 2389: return 191.962600221;    // 192Ir
		case 2390: return 192.962921621;    // 193Ir
		case 2391: return 193.965073521;    // 194Ir
		case 2392: return 194.965974721;    // 195Ir
		case 2393: return 195.96839741;     // 196Ir
		case 2394: return 196.96965522;     // 197Ir
		case 2395: return 197.9722821;      // 198Ir
		case 2396: return 198.97380544;     // 199Ir
		case 2397: return 199.9768021;      // 200Ir
		case 2398: return 200.9786421;      // 201Ir
		case 2399: return 201.9819932;      // 202Ir
		case 2400: return 202.9842343;      // 203Ir
		case 2401: return 203.9896043;      // 204Ir
		case 2402: return 165.9948654;      // 166Pt
		case 2403: return 166.9926933;      // 167Pt
		case 2404: return 167.9881317;      // 168Pt
		case 2405: return 168.9865722;      // 169Pt
		case 2406: return 169.98249620;     // 170Pt
		case 2407: return 170.98124578;     // 171Pt
		case 2408: return 171.97735112;     // 172Pt
		case 2409: return 172.97644360;     // 173Pt
		case 2410: return 173.97282011;     // 174Pt
		case 2411: return 174.97241019;     // 175Pt
		case 2412: return 175.96893814;     // 176Pt
		case 2413: return 176.96847016;     // 177Pt
		case 2414: return 177.96565011;     // 178Pt
		case 2415: return 178.965359086;    // 179Pt
		case 2416: return 179.96303212;     // 180Pt
		case 2417: return 180.96309816;     // 181Pt
		case 2418: return 181.96117214;     // 182Pt
		case 2419: return 182.96159717;     // 183Pt
		case 2420: return 183.95991517;     // 184Pt
		case 2421: return 184.96061428;     // 185Pt
		case 2422: return 185.95935123;     // 186Pt
		case 2423: return 186.96061726;     // 187Pt
		case 2424: return 187.959388961;    // 188Pt
		case 2425: return 188.96083112;     // 189Pt
		case 2426: return 189.959929763;    // 190Pt
		case 2427: return 190.961672953;    // 191Pt
		case 2428: return 191.961038732;    // 192Pt
		case 2429: return 192.962982421;    // 193Pt
		case 2430: return 193.962680910;    // 194Pt
		case 2431: return 194.964791710;    // 195Pt
		case 2432: return 195.9649520999;   // 196Pt
		case 2433: return 196.9673406994;   // 197Pt
		case 2434: return 197.967894923;    // 198Pt
		case 2435: return 198.970595224;    // 199Pt
		case 2436: return 199.97144322;     // 200Pt
		case 2437: return 200.97451354;     // 201Pt
		case 2438: return 201.97563927;     // 202Pt
		case 2439: return 202.9789321;      // 203Pt
		case 2440: return 203.9807621;      // 204Pt
		case 2441: return 204.9860832;      // 205Pt
		case 2442: return 205.9896632;      // 206Pt
		case 2443: return 168.9980832;      // 169Au
		case 2444: return 169.9959722;      // 170Au
		case 2445: return 170.99187622;     // 171Au
		case 2446: return 171.98994281;     // 172Au
		case 2447: return 172.98624126;     // 173Au
		case 2448: return 173.98471795;     // 174Au
		case 2449: return 174.98130442;     // 175Au
		case 2450: return 175.98025036;     // 176Au
		case 2451: return 176.97687011;     // 177Au
		case 2452: return 177.97603261;     // 178Au
		case 2453: return 178.97317413;     // 179Au
		case 2454: return 179.97252321;     // 180Au
		case 2455: return 180.97007921;     // 181Au
		case 2456: return 181.96961822;     // 182Au
		case 2457: return 182.96759110;     // 183Au
		case 2458: return 183.96745224;     // 184Au
		case 2459: return 184.96579028;     // 185Au
		case 2460: return 185.96595323;     // 186Au
		case 2461: return 186.96454324;     // 187Au
		case 2462: return 187.96534917;     // 188Au
		case 2463: return 188.96394822;     // 189Au
		case 2464: return 189.96469817;     // 190Au
		case 2465: return 190.96370240;     // 191Au
		case 2466: return 191.96481417;     // 192Au
		case 2467: return 192.964137393;    // 193Au
		case 2468: return 193.965417823;    // 194Au
		case 2469: return 194.965035215;    // 195Au
		case 2470: return 195.966569932;    // 196Au
		case 2471: return 196.9665687971;   // 197Au
		case 2472: return 197.9682424270;   // 198Au
		case 2473: return 198.9687652870;   // 199Au
		case 2474: return 199.97075629;     // 200Au
		case 2475: return 200.971657534;    // 201Au
		case 2476: return 201.97385625;     // 202Au
		case 2477: return 202.975154433;    // 203Au
		case 2478: return 203.9778322;      // 204Au
		case 2479: return 204.9798521;      // 205Au
		case 2480: return 205.9847432;      // 206Au
		case 2481: return 206.9884032;      // 207Au
		case 2482: return 207.9934532;      // 208Au
		case 2483: return 208.9973543;      // 209Au
		case 2484: return 210.0025043;      // 210Au
		case 2485: return 171.0035333;      // 171Hg
		case 2486: return 171.9988117;      // 172Hg
		case 2487: return 172.9970922;      // 173Hg
		case 2488: return 173.99286521;     // 174Hg
		case 2489: return 174.99144178;     // 175Hg
		case 2490: return 175.98736114;     // 176Hg
		case 2491: return 176.98627781;     // 177Hg
		case 2492: return 177.98248412;     // 178Hg
		case 2493: return 178.98183129;     // 179Hg
		case 2494: return 179.97826014;     // 180Hg
		case 2495: return 180.97781917;     // 181Hg
		case 2496: return 181.97468911;     // 182Hg
		case 2497: return 182.974444876;    // 183Hg
		case 2498: return 183.97171411;     // 184Hg
		case 2499: return 184.97189917;     // 185Hg
		case 2500: return 185.96936213;     // 186Hg
		case 2501: return 186.96981415;     // 187Hg
		case 2502: return 187.96756712;     // 188Hg
		case 2503: return 188.96819534;     // 189Hg
		case 2504: return 189.96632317;     // 190Hg
		case 2505: return 190.96715724;     // 191Hg
		case 2506: return 191.96563517;     // 192Hg
		case 2507: return 192.96665317;     // 193Hg
		case 2508: return 193.965449131;    // 194Hg
		case 2509: return 194.96672125;     // 195Hg
		case 2510: return 195.965832632;    // 196Hg
		case 2511: return 196.967212835;    // 197Hg
		case 2512: return 197.9667686052;   // 198Hg
		case 2513: return 198.9682806446;   // 199Hg
		case 2514: return 199.9683265947;   // 200Hg
		case 2515: return 200.9703028469;   // 201Hg
		case 2516: return 201.9706434069;   // 202Hg
		case 2517: return 202.972872818;    // 203Hg
		case 2518: return 203.9734939853;   // 204Hg
		case 2519: return 204.976073439;    // 205Hg
		case 2520: return 205.97751422;     // 206Hg
		case 2521: return 206.98230032;     // 207Hg
		case 2522: return 207.98575933;     // 208Hg
		case 2523: return 208.9907216;      // 209Hg
		case 2524: return 209.9942421;      // 210Hg
		case 2525: return 210.9993321;      // 211Hg
		case 2526: return 212.0029632;      // 212Hg
		case 2527: return 213.0082332;      // 213Hg
		case 2528: return 214.0120043;      // 214Hg
		case 2529: return 215.0174043;      // 215Hg
		case 2530: return 216.0213243;      // 216Hg
		case 2531: return 176.00062481;     // 176Tl
		case 2532: return 176.99643125;     // 177Tl
		case 2533: return 177.9948511;      // 178Tl
		case 2534: return 178.99111143;     // 179Tl
		case 2535: return 179.99005764;     // 180Tl
		case 2536: return 180.986260098;    // 181Tl
		case 2537: return 181.98571363;     // 182Tl
		case 2538: return 182.98219310;     // 183Tl
		case 2539: return 183.98188622;     // 184Tl
		case 2540: return 184.97878922;     // 185Tl
		case 2541: return 185.97865124;     // 186Tl
		case 2542: return 186.975906388;    // 187Tl
		case 2543: return 187.97602132;     // 188Tl
		case 2544: return 188.97358812;     // 189Tl
		case 2545: return 189.97382854;     // 190Tl
		case 2546: return 190.971784279;    // 191Tl
		case 2547: return 191.97222534;     // 192Tl
		case 2548: return 192.970502072;    // 193Tl
		case 2549: return 193.97108115;     // 194Tl
		case 2550: return 194.96977412;     // 195Tl
		case 2551: return 195.97048113;     // 196Tl
		case 2552: return 196.96957618;     // 197Tl
		case 2553: return 197.97048386;     // 198Tl
		case 2554: return 198.96987730;     // 199Tl
		case 2555: return 199.970963362;    // 200Tl
		case 2556: return 200.97082215;     // 201Tl
		case 2557: return 201.97210215;     // 202Tl
		case 2558: return 202.972344614;    // 203Tl
		case 2559: return 203.973863913;    // 204Tl
		case 2560: return 204.974427814;    // 205Tl
		case 2561: return 205.976110615;    // 206Tl
		case 2562: return 206.977419759;    // 207Tl
		case 2563: return 207.982019021;    // 208Tl
		case 2564: return 208.985359486;    // 209Tl
		case 2565: return 209.99007412;     // 210Tl
		case 2566: return 210.99347545;     // 211Tl
		case 2567: return 211.9983422;      // 212Tl
		case 2568: return 213.00191529;     // 213Tl
		case 2569: return 214.0069421;      // 214Tl
		case 2570: return 215.0106432;      // 215Tl
		case 2571: return 216.0158032;      // 216Tl
		case 2572: return 217.0196643;      // 217Tl
		case 2573: return 218.0247943;      // 218Tl
		case 2574: return 178.00383126;     // 178Pb
		case 2575: return 179.00220181;     // 179Pb
		case 2576: return 179.99792815;     // 180Pb
		case 2577: return 180.99665381;     // 181Pb
		case 2578: return 181.99267213;     // 182Pb
		case 2579: return 182.99187230;     // 183Pb
		case 2580: return 183.98813614;     // 184Pb
		case 2581: return 184.98761017;     // 185Pb
		case 2582: return 185.98423812;     // 186Pb
		case 2583: return 186.983910955;    // 187Pb
		case 2584: return 187.98087511;     // 188Pb
		case 2585: return 188.98080737;     // 189Pb
		case 2586: return 189.97808213;     // 190Pb
		case 2587: return 190.97827641;     // 191Pb
		case 2588: return 191.97577513;     // 192Pb
		case 2589: return 192.97617353;     // 193Pb
		case 2590: return 193.97401219;     // 194Pb
		case 2591: return 194.97454325;     // 195Pb
		case 2592: return 195.97277415;     // 196Pb
		case 2593: return 196.973431260;    // 197Pb
		case 2594: return 197.97203416;     // 198Pb
		case 2595: return 198.97291311;     // 199Pb
		case 2596: return 199.97181912;     // 200Pb
		case 2597: return 200.97288323;     // 201Pb
		case 2598: return 201.972152040;    // 202Pb
		case 2599: return 202.973391171;    // 203Pb
		case 2600: return 203.973044013;    // 204Pb
		case 2601: return 204.974482213;    // 205Pb
		case 2602: return 205.974465713;    // 206Pb
		case 2603: return 206.975897313;    // 207Pb
		case 2604: return 207.976652513;    // 208Pb
		case 2605: return 208.981090519;    // 209Pb
		case 2606: return 209.984188916;    // 210Pb
		case 2607: return 210.988737128;    // 211Pb
		case 2608: return 211.991897723;    // 212Pb
		case 2609: return 212.996562972;    // 213Pb
		case 2610: return 213.999805925;    // 214Pb
		case 2611: return 215.0047411;      // 215Pb
		case 2612: return 216.0080321;      // 216Pb
		case 2613: return 217.0131432;      // 217Pb
		case 2614: return 218.0165932;      // 218Pb
		case 2615: return 219.0217743;      // 219Pb
		case 2616: return 220.0254143;      // 220Pb
		case 2617: return 184.00127584;     // 184Bi
		case 2618: return 184.99760087;     // 185Bi
		case 2619: return 185.99664465;     // 186Bi
		case 2620: return 186.99314711;     // 187Bi
		case 2621: return 187.99228722;     // 188Bi
		case 2622: return 188.98919522;     // 189Bi
		case 2623: return 189.98862224;     // 190Bi
		case 2624: return 190.985786680;    // 191Bi
		case 2625: return 191.98546933;     // 192Bi
		case 2626: return 192.98296010;     // 193Bi
		case 2627: return 193.98278554;     // 194Bi
		case 2628: return 194.980648857;    // 195Bi
		case 2629: return 195.98066726;     // 196Bi
		case 2630: return 196.978865189;    // 197Bi
		case 2631: return 197.97920630;     // 198Bi
		case 2632: return 198.97767311;     // 199Bi
		case 2633: return 199.97813124;     // 200Bi
		case 2634: return 200.97701016;     // 201Bi
		case 2635: return 201.97773417;     // 202Bi
		case 2636: return 202.97689314;     // 203Bi
		case 2637: return 203.977836199;    // 204Bi
		case 2638: return 204.977386755;    // 205Bi
		case 2639: return 205.978499382;    // 206Bi
		case 2640: return 206.978471026;    // 207Bi
		case 2641: return 207.979742525;    // 208Bi
		case 2642: return 208.980399116;    // 209Bi
		case 2643: return 209.984120716;    // 210Bi
		case 2644: return 210.987269759;    // 211Bi
		case 2645: return 211.991286021;    // 212Bi
		case 2646: return 212.994385156;    // 213Bi
		case 2647: return 213.99871212;     // 214Bi
		case 2648: return 215.00177016;     // 215Bi
		case 2649: return 216.00630612;     // 216Bi
		case 2650: return 217.00937219;     // 217Bi
		case 2651: return 218.01418829;     // 218Bi
		case 2652: return 219.0174821;      // 219Bi
		case 2653: return 220.0223532;      // 220Bi
		case 2654: return 221.0258732;      // 221Bi
		case 2655: return 222.0307832;      // 222Bi
		case 2656: return 223.0345043;      // 223Bi
		case 2657: return 224.0394743;      // 224Bi
		case 2658: return 186.00439335;     // 186Po
		case 2659: return 187.00304134;     // 187Po
		case 2660: return 187.99941621;     // 188Po
		case 2661: return 188.99847324;     // 189Po
		case 2662: return 189.99510114;     // 190Po
		case 2663: return 190.994558576;    // 191Po
		case 2664: return 191.99133612;     // 192Po
		case 2665: return 192.99102637;     // 193Po
		case 2666: return 193.98818614;     // 194Po
		case 2667: return 194.98812641;     // 195Po
		case 2668: return 195.98552614;     // 196Po
		case 2669: return 196.98566053;     // 197Po
		case 2670: return 197.98338919;     // 198Po
		case 2671: return 198.98366725;     // 199Po
		case 2672: return 199.98179915;     // 200Po
		case 2673: return 200.982259863;    // 201Po
		case 2674: return 201.98075816;     // 202Po
		case 2675: return 202.981416193;    // 203Po
		case 2676: return 203.98031012;     // 204Po
		case 2677: return 204.98120322;     // 205Po
		case 2678: return 205.980474043;    // 206Po
		case 2679: return 206.981593872;    // 207Po
		case 2680: return 207.981246119;    // 208Po
		case 2681: return 208.982430820;    // 209Po
		case 2682: return 209.982874113;    // 210Po
		case 2683: return 210.986653614;    // 211Po
		case 2684: return 211.988868413;    // 212Po
		case 2685: return 212.992857633;    // 213Po
		case 2686: return 213.995201716;    // 214Po
		case 2687: return 214.999420127;    // 215Po
		case 2688: return 216.001915223;    // 216Po
		case 2689: return 217.006318267;    // 217Po
		case 2690: return 218.008973525;    // 218Po
		case 2691: return 219.01361417;     // 219Po
		case 2692: return 220.01638619;     // 220Po
		case 2693: return 221.02122821;     // 221Po
		case 2694: return 222.02414043;     // 222Po
		case 2695: return 223.0290721;      // 223Po
		case 2696: return 224.0321121;      // 224Po
		case 2697: return 225.0370732;      // 225Po
		case 2698: return 226.0403143;      // 226Po
		case 2699: return 227.0453943;      // 227Po
		case 2700: return 191.00414817;     // 191At
		case 2701: return 192.00315235;     // 192At
		case 2702: return 192.99992723;     // 193At
		case 2703: return 193.99923629;     // 194At
		case 2704: return 194.996268598;    // 195At
		case 2705: return 195.99580033;     // 196At
		case 2706: return 196.99318955;     // 197At
		case 2707: return 197.99278454;     // 198At
		case 2708: return 198.990527758;    // 199At
		case 2709: return 199.99035126;     // 200At
		case 2710: return 200.988417188;    // 201At
		case 2711: return 201.98863030;     // 202At
		case 2712: return 202.98694311;     // 203At
		case 2713: return 203.98725124;     // 204At
		case 2714: return 204.98607616;     // 205At
		case 2715: return 205.98665716;     // 206At
		case 2716: return 206.98580013;     // 207At
		case 2717: return 207.986613396;    // 208At
		case 2718: return 208.986170255;    // 209At
		case 2719: return 209.987147983;    // 210At
		case 2720: return 210.987496630;    // 211At
		case 2721: return 211.990737726;    // 212At
		case 2722: return 212.992937053;    // 213At
		case 2723: return 213.996372146;    // 214At
		case 2724: return 214.998652873;    // 215At
		case 2725: return 216.002423639;    // 216At
		case 2726: return 217.004719255;    // 217At
		case 2727: return 218.00869512;     // 218At
		case 2728: return 219.011161842;    // 219At
		case 2729: return 220.01543315;     // 220At
		case 2730: return 221.01801715;     // 221At
		case 2731: return 222.02249417;     // 222At
		case 2732: return 223.02515115;     // 223At
		case 2733: return 224.02974924;     // 224At
		case 2734: return 225.0326332;      // 225At
		case 2735: return 226.0371632;      // 226At
		case 2736: return 227.0402432;      // 227At
		case 2737: return 228.0447543;      // 228At
		case 2738: return 229.0481243;      // 229At
		case 2739: return 193.00970827;     // 193Rn
		case 2740: return 194.00614418;     // 194Rn
		case 2741: return 195.00542254;     // 195Rn
		case 2742: return 196.00211615;     // 196Rn
		case 2743: return 197.00158538;     // 197Rn
		case 2744: return 197.99867914;     // 198Rn
		case 2745: return 198.99839068;     // 199Rn
		case 2746: return 199.99569014;     // 200Rn
		case 2747: return 200.99562853;     // 201Rn
		case 2748: return 201.99326419;     // 202Rn
		case 2749: return 202.99338825;     // 203Rn
		case 2750: return 203.99143016;     // 204Rn
		case 2751: return 204.99171954;     // 205Rn
		case 2752: return 205.99021416;     // 206Rn
		case 2753: return 206.990730391;    // 207Rn
		case 2754: return 207.98963512;     // 208Rn
		case 2755: return 208.99041522;     // 209Rn
		case 2756: return 209.989689149;    // 210Rn
		case 2757: return 210.990601173;    // 211Rn
		case 2758: return 211.990703934;    // 212Rn
		case 2759: return 212.993883161;    // 213Rn
		case 2760: return 213.995363099;    // 214Rn
		case 2761: return 214.998745983;    // 215Rn
		case 2762: return 216.000271965;    // 216Rn
		case 2763: return 217.003928045;    // 217Rn
		case 2764: return 218.005601625;    // 218Rn
		case 2765: return 219.009480427;    // 219Rn
		case 2766: return 220.011394123;    // 220Rn
		case 2767: return 221.015537163;    // 221Rn
		case 2768: return 222.017578225;    // 222Rn
		case 2769: return 223.021889384;    // 223Rn
		case 2770: return 224.02409611;     // 224Rn
		case 2771: return 225.02848612;     // 225Rn
		case 2772: return 226.03086111;     // 226Rn
		case 2773: return 227.03530415;     // 227Rn
		case 2774: return 228.03783519;     // 228Rn
		case 2775: return 229.04225714;     // 229Rn
		case 2776: return 230.0451421;      // 230Rn
		case 2777: return 231.0498732;      // 231Rn
		case 2778: return 199.00725945;     // 199Fr
		case 2779: return 200.00658663;     // 200Fr
		case 2780: return 201.00386777;     // 201Fr
		case 2781: return 202.00332055;     // 202Fr
		case 2782: return 203.000940767;    // 203Fr
		case 2783: return 204.00065226;     // 204Fr
		case 2784: return 204.998593984;    // 205Fr
		case 2785: return 205.99866630;     // 206Fr
		case 2786: return 206.99694619;     // 207Fr
		case 2787: return 207.99713812;     // 208Fr
		case 2788: return 208.99595516;     // 209Fr
		case 2789: return 209.99642216;     // 210Fr
		case 2790: return 210.99555613;     // 211Fr
		case 2791: return 211.996225794;    // 212Fr
		case 2792: return 212.996186055;    // 213Fr
		case 2793: return 213.998971393;    // 214Fr
		case 2794: return 215.000341876;    // 215Fr
		case 2795: return 216.003189945;    // 216Fr
		case 2796: return 217.004632370;    // 217Fr
		case 2797: return 218.007578751;    // 218Fr
		case 2798: return 219.009252476;    // 219Fr
		case 2799: return 220.012327744;    // 220Fr
		case 2800: return 221.014255254;    // 221Fr
		case 2801: return 222.01755223;     // 222Fr
		case 2802: return 223.019736025;    // 223Fr
		case 2803: return 224.02339814;     // 224Fr
		case 2804: return 225.02557313;     // 225Fr
		case 2805: return 226.02956613;     // 226Fr
		case 2806: return 227.03186914;     // 227Fr
		case 2807: return 228.03582314;     // 228Fr
		case 2808: return 229.03829815;     // 229Fr
		case 2809: return 230.04241617;     // 230Fr
		case 2810: return 231.04515827;     // 231Fr
		case 2811: return 232.0493717;      // 232Fr
		case 2812: return 233.0526432;      // 233Fr
		case 2813: return 201.0127111;      // 201Ra
		case 2814: return 202.00976026;     // 202Ra
		case 2815: return 203.00930486;     // 203Ra
		case 2816: return 204.00649216;     // 204Ra
		case 2817: return 205.00626876;     // 205Ra
		case 2818: return 206.00382819;     // 206Ra
		case 2819: return 207.00379959;     // 207Ra
		case 2820: return 208.00184117;     // 208Ra
		case 2821: return 209.00199054;     // 209Ra
		case 2822: return 210.00049416;     // 210Ra
		case 2823: return 211.000893285;    // 211Ra
		case 2824: return 211.99978712;     // 212Ra
		case 2825: return 213.00038422;     // 213Ra
		case 2826: return 214.000099756;    // 214Ra
		case 2827: return 215.002720482;    // 215Ra
		case 2828: return 216.003533494;    // 216Ra
		case 2829: return 217.006320792;    // 217Ra
		case 2830: return 218.00714112;     // 218Ra
		case 2831: return 219.010085589;    // 219Ra
		case 2832: return 220.011025989;    // 220Ra
		case 2833: return 221.013917750;    // 221Ra
		case 2834: return 222.015374849;    // 222Ra
		case 2835: return 223.018502327;    // 223Ra
		case 2836: return 224.020212023;    // 224Ra
		case 2837: return 225.023611932;    // 225Ra
		case 2838: return 226.025410325;    // 226Ra
		case 2839: return 227.029178325;    // 227Ra
		case 2840: return 228.031070726;    // 228Ra
		case 2841: return 229.03494216;     // 229Ra
		case 2842: return 230.03705511;     // 230Ra
		case 2843: return 231.04102712;     // 231Ra
		case 2844: return 232.043475398;    // 232Ra
		case 2845: return 233.04758217;     // 233Ra
		case 2846: return 234.05034233;     // 234Ra
		case 2847: return 235.0549732;      // 235Ra
		case 2848: return 206.01445277;     // 206Ac
		case 2849: return 207.01196654;     // 207Ac
		case 2850: return 208.01155060;     // 208Ac
		case 2851: return 209.00949554;     // 209Ac
		case 2852: return 210.00943662;     // 210Ac
		case 2853: return 211.00773257;     // 211Ac
		case 2854: return 212.00781355;     // 212Ac
		case 2855: return 213.00660956;     // 213Ac
		case 2856: return 214.00691816;     // 214Ac
		case 2857: return 215.00647513;     // 215Ac
		case 2858: return 216.00874312;     // 216Ac
		case 2859: return 217.00934412;     // 217Ac
		case 2860: return 218.01164254;     // 218Ac
		case 2861: return 219.01242154;     // 219Ac
		case 2862: return 220.014754966;    // 220Ac
		case 2863: return 221.01559254;     // 221Ac
		case 2864: return 222.017844256;    // 222Ac
		case 2865: return 223.019137777;    // 223Ac
		case 2866: return 224.021723245;    // 224Ac
		case 2867: return 225.023230053;    // 225Ac
		case 2868: return 226.026098436;    // 226Ac
		case 2869: return 227.027752325;    // 227Ac
		case 2870: return 228.031021527;    // 228Ac
		case 2871: return 229.03295613;     // 229Ac
		case 2872: return 230.03632717;     // 230Ac
		case 2873: return 231.03839314;     // 231Ac
		case 2874: return 232.04203414;     // 232Ac
		case 2875: return 233.04434614;     // 233Ac
		case 2876: return 234.04813915;     // 234Ac
		case 2877: return 235.05084015;     // 235Ac
		case 2878: return 236.05498841;     // 236Ac
		case 2879: return 237.0582743;      // 237Ac
		case 2880: return 208.01790036;     // 208Th
		case 2881: return 209.01775393;     // 209Th
		case 2882: return 210.01509420;     // 210Th
		case 2883: return 211.01492980;     // 211Th
		case 2884: return 212.01298817;     // 212Th
		case 2885: return 213.01300976;     // 213Th
		case 2886: return 214.01150017;     // 214Th
		case 2887: return 215.011724895;    // 215Th
		case 2888: return 216.01105613;     // 216Th
		case 2889: return 217.01311722;     // 217Th
		case 2890: return 218.01327611;     // 218Th
		case 2891: return 219.01553754;     // 219Th
		case 2892: return 220.01574824;     // 220Th
		case 2893: return 221.01818410;     // 221Th
		case 2894: return 222.01846913;     // 222Th
		case 2895: return 223.020811999;    // 223Th
		case 2896: return 224.02146411;     // 224Th
		case 2897: return 225.023951455;    // 225Th
		case 2898: return 226.024903450;    // 226Th
		case 2899: return 227.027704227;    // 227Th
		case 2900: return 228.028741323;    // 228Th
		case 2901: return 229.031762730;    // 229Th
		case 2902: return 230.033134119;    // 230Th
		case 2903: return 231.036304619;    // 231Th
		case 2904: return 232.038055821;    // 232Th
		case 2905: return 233.041582321;    // 233Th
		case 2906: return 234.043601437;    // 234Th
		case 2907: return 235.04725514;     // 235Th
		case 2908: return 236.04965715;     // 236Th
		case 2909: return 237.05362917;     // 237Th
		case 2910: return 238.0565030;      // 238Th
		case 2911: return 239.0607743;      // 239Th
		case 2912: return 212.02320380;     // 212Pa
		case 2913: return 213.02110976;     // 213Pa
		case 2914: return 214.02091882;     // 214Pa
		case 2915: return 215.01918378;     // 215Pa
		case 2916: return 216.01910957;     // 216Pa
		case 2917: return 217.01832556;     // 217Pa
		case 2918: return 218.02005920;     // 218Pa
		case 2919: return 219.01990455;     // 219Pa
		case 2920: return 220.02170555;     // 220Pa
		case 2921: return 221.02187555;     // 221Pa
		case 2922: return 222.02378478;     // 222Pa
		case 2923: return 223.02396376;     // 223Pa
		case 2924: return 224.025617682;    // 224Pa
		case 2925: return 225.02613176;     // 225Pa
		case 2926: return 226.02794812;     // 226Pa
		case 2927: return 227.028805480;    // 227Pa
		case 2928: return 228.031051747;    // 228Pa
		case 2929: return 229.032097238;    // 229Pa
		case 2930: return 230.034541035;    // 230Pa
		case 2931: return 231.035884224;    // 231Pa
		case 2932: return 232.038591783;    // 232Pa
		case 2933: return 233.040247222;    // 233Pa
		case 2934: return 234.043307251;    // 234Pa
		case 2935: return 235.04539915;     // 235Pa
		case 2936: return 236.04866815;     // 236Pa
		case 2937: return 237.05102314;     // 237Pa
		case 2938: return 238.05463717;     // 238Pa
		case 2939: return 239.0572621;      // 239Pa
		case 2940: return 240.0609832;      // 240Pa
		case 2941: return 241.0640843;      // 241Pa
		case 2942: return 217.0246611;      // 217U
		case 2943: return 218.02352320;     // 218U
		case 2944: return 219.02499955;     // 219U
		case 2945: return 220.0246211;      // 220U
		case 2946: return 221.0262811;      // 221U
		case 2947: return 222.0260011;      // 222U
		case 2948: return 223.02773976;     // 223U
		case 2949: return 224.02760527;     // 224U
		case 2950: return 225.02939113;     // 225U
		case 2951: return 226.02933914;     // 226U
		case 2952: return 227.03115718;     // 227U
		case 2953: return 228.03137115;     // 228U
		case 2954: return 229.033506364;    // 229U
		case 2955: return 230.033940151;    // 230U
		case 2956: return 231.036293932;    // 231U
		case 2957: return 232.037156323;    // 232U
		case 2958: return 233.039635529;    // 233U
		case 2959: return 234.040952319;    // 234U
		case 2960: return 235.043930119;    // 235U
		case 2961: return 236.045568219;    // 236U
		case 2962: return 237.048730420;    // 237U
		case 2963: return 238.050788420;    // 238U
		case 2964: return 239.054293520;    // 239U
		case 2965: return 240.056593457;    // 240U
		case 2966: return 241.0603332;      // 241U
		case 2967: return 242.0629322;      // 242U
		case 2968: return 243.0669943;      // 243U
		case 2969: return 219.0314321;      // 219Np
		case 2970: return 220.0325421;      // 220Np
		case 2971: return 221.0320421;      // 221Np
		case 2972: return 222.0333021;      // 222Np
		case 2973: return 223.0328521;      // 223Np
		case 2974: return 224.0342221;      // 224Np
		case 2975: return 225.03391177;     // 225Np
		case 2976: return 226.03518895;     // 226Np
		case 2977: return 227.03495778;     // 227Np
		case 2978: return 228.03606754;     // 228Np
		case 2979: return 229.03626493;     // 229Np
		case 2980: return 230.03782855;     // 230Np
		case 2981: return 231.03824554;     // 231Np
		case 2982: return 232.0401111;      // 232Np
		case 2983: return 233.04074155;     // 233Np
		case 2984: return 234.042895391;    // 234Np
		case 2985: return 235.044063521;    // 235Np
		case 2986: return 236.04657054;     // 236Np
		case 2987: return 237.048173619;    // 237Np
		case 2988: return 238.050946619;    // 238Np
		case 2989: return 239.052939222;    // 239Np
		case 2990: return 240.05616518;     // 240Np
		case 2991: return 241.05825376;     // 241Np
		case 2992: return 242.0616421;      // 242Np
		case 2993: return 243.06428034;     // 243Np
		case 2994: return 244.0678532;      // 244Np
		case 2995: return 245.0708043;      // 245Np
		case 2996: return 228.03873233;     // 228Pu
		case 2997: return 229.04014455;     // 229Pu
		case 2998: return 230.03965016;     // 230Pu
		case 2999: return 231.04110228;     // 231Pu
		case 3000: return 232.04118519;     // 232Pu
		case 3001: return 233.04299854;     // 233Pu
		case 3002: return 234.043317475;    // 234Pu
		case 3003: return 235.04528622;     // 235Pu
		case 3004: return 236.046058123;    // 236Pu
		case 3005: return 237.048409824;    // 237Pu
		case 3006: return 238.049560119;    // 238Pu
		case 3007: return 239.052163619;    // 239Pu
		case 3008: return 240.053813819;    // 240Pu
		case 3009: return 241.056851719;    // 241Pu
		case 3010: return 242.058742820;    // 242Pu
		case 3011: return 243.062003634;    // 243Pu
		case 3012: return 244.064205356;    // 244Pu
		case 3013: return 245.06782615;     // 245Pu
		case 3014: return 246.07020516;     // 246Pu
		case 3015: return 247.0741921;      // 247Pu
		case 3016: return 230.0460914;      // 230Am
		case 3017: return 231.0455632;      // 231Am
		case 3018: return 232.0464532;      // 232Am
		case 3019: return 233.0464411;      // 233Am
		case 3020: return 234.0477317;      // 234Am
		case 3021: return 235.04790856;     // 235Am
		case 3022: return 236.0494312;      // 236Am
		case 3023: return 237.04999664;     // 237Am
		case 3024: return 238.05198554;     // 238Am
		case 3025: return 239.053024726;    // 239Am
		case 3026: return 240.05530015;     // 240Am
		case 3027: return 241.056829319;    // 241Am
		case 3028: return 242.059549419;    // 242Am
		case 3029: return 243.061381324;    // 243Am
		case 3030: return 244.064285122;    // 244Am
		case 3031: return 245.066454834;    // 245Am
		case 3032: return 246.06977520;     // 246Am
		case 3033: return 247.0720911;      // 247Am
		case 3034: return 248.0757522;      // 248Am
		case 3035: return 249.0784832;      // 249Am
		case 3036: return 232.0498222;      // 232Cm
		case 3037: return 233.05077077;     // 233Cm
		case 3038: return 234.05016020;     // 234Cm
		case 3039: return 235.0515422;      // 235Cm
		case 3040: return 236.05137420;     // 236Cm
		case 3041: return 237.05286976;     // 237Cm
		case 3042: return 238.05308113;     // 238Cm
		case 3043: return 239.05491058;     // 239Cm
		case 3044: return 240.055529724;    // 240Cm
		case 3045: return 241.057653223;    // 241Cm
		case 3046: return 242.058836019;    // 242Cm
		case 3047: return 243.061389322;    // 243Cm
		case 3048: return 244.062752819;    // 244Cm
		case 3049: return 245.065491522;    // 245Cm
		case 3050: return 246.067223822;    // 246Cm
		case 3051: return 247.070354147;    // 247Cm
		case 3052: return 248.072349956;    // 248Cm
		case 3053: return 249.075954856;    // 249Cm
		case 3054: return 250.07835812;     // 250Cm
		case 3055: return 251.08228624;     // 251Cm
		case 3056: return 252.0848732;      // 252Cm
		case 3057: return 234.0572715;      // 234Bk
		case 3058: return 235.0565843;      // 235Bk
		case 3059: return 236.0574843;      // 236Bk
		case 3060: return 237.0571024;      // 237Bk
		case 3061: return 238.0582027;      // 238Bk
		case 3062: return 239.0582422;      // 239Bk
		case 3063: return 240.0597616;      // 240Bk
		case 3064: return 241.0601622;      // 241Bk
		case 3065: return 242.0619822;      // 242Bk
		case 3066: return 243.063007851;    // 243Bk
		case 3067: return 244.06518116;     // 244Bk
		case 3068: return 245.066361824;    // 245Bk
		case 3069: return 246.06867364;     // 246Bk
		case 3070: return 247.070307359;    // 247Bk
		case 3071: return 248.07308876;     // 248Bk
		case 3072: return 249.074987727;    // 249Bk
		case 3073: return 250.078316742;    // 250Bk
		case 3074: return 251.08076212;     // 251Bk
		case 3075: return 252.0843122;      // 252Bk
		case 3076: return 253.0868839;      // 253Bk
		case 3077: return 254.0906032;      // 254Bk
		case 3078: return 237.06219894;     // 237Cf
		case 3079: return 238.0614932;      // 238Cf
		case 3080: return 239.0625323;      // 239Cf
		case 3081: return 240.06225620;     // 240Cf
		case 3082: return 241.0636918;      // 241Cf
		case 3083: return 242.06375414;     // 242Cf
		case 3084: return 243.0654812;      // 243Cf
		case 3085: return 244.066000831;    // 244Cf
		case 3086: return 245.068048730;    // 245Cf
		case 3087: return 246.068805522;    // 246Cf
		case 3088: return 247.07096516;     // 247Cf
		case 3089: return 248.072185157;    // 248Cf
		case 3090: return 249.074853923;    // 249Cf
		case 3091: return 250.076406222;    // 250Cf
		case 3092: return 251.079588648;    // 251Cf
		case 3093: return 252.081627256;    // 252Cf
		case 3094: return 253.085134567;    // 253Cf
		case 3095: return 254.08732413;     // 254Cf
		case 3096: return 255.0910522;      // 255Cf
		case 3097: return 256.0934434;      // 256Cf
		case 3098: return 239.0682332;      // 239Es
		case 3099: return 240.0689243;      // 240Es
		case 3100: return 241.0685624;      // 241Es
		case 3101: return 242.0695728;      // 242Es
		case 3102: return 243.0695122;      // 243Es
		case 3103: return 244.0708820;      // 244Es
		case 3104: return 245.0712522;      // 245Es
		case 3105: return 246.0729024;      // 246Es
		case 3106: return 247.07362221;     // 247Es
		case 3107: return 248.07547156;     // 248Es
		case 3108: return 249.07641132;     // 249Es
		case 3109: return 250.0786111;      // 250Es
		case 3110: return 251.079993667;    // 251Es
		case 3111: return 252.08298054;     // 252Es
		case 3112: return 253.084825727;    // 253Es
		case 3113: return 254.088022245;    // 254Es
		case 3114: return 255.09027512;     // 255Es
		case 3115: return 256.0936011;      // 256Es
		case 3116: return 257.0959844;      // 257Es
		case 3117: return 258.0995232;      // 258Es
		case 3118: return 241.0742132;      // 241Fm
		case 3119: return 242.0734343;      // 242Fm
		case 3120: return 243.0744623;      // 243Fm
		case 3121: return 244.0740422;      // 244Fm
		case 3122: return 245.0753521;      // 245Fm
		case 3123: return 246.07535017;     // 246Fm
		case 3124: return 247.0769412;      // 247Fm
		case 3125: return 248.077186592;    // 248Fm
		case 3126: return 249.078927568;    // 249Fm
		case 3127: return 250.079521086;    // 250Fm
		case 3128: return 251.08154016;     // 251Fm
		case 3129: return 252.082467161;    // 252Fm
		case 3130: return 253.085184637;    // 253Fm
		case 3131: return 254.086854430;    // 254Fm
		case 3132: return 255.089964052;    // 255Fm
		case 3133: return 256.091774578;    // 256Fm
		case 3134: return 257.095106169;    // 257Fm
		case 3135: return 258.0970822;      // 258Fm
		case 3136: return 259.1006030;      // 259Fm
		case 3137: return 260.1028155;      // 260Fm
		case 3138: return 245.0808133;      // 245Md
		case 3139: return 246.0817128;      // 246Md
		case 3140: return 247.0815222;      // 247Md
		case 3141: return 248.0828226;      // 248Md
		case 3142: return 249.0829122;      // 249Md
		case 3143: return 250.0844132;      // 250Md
		case 3144: return 251.08477420;     // 251Md
		case 3145: return 252.0864314;      // 252Md
		case 3146: return 253.08714434;     // 253Md
		case 3147: return 254.0895911;      // 254Md
		case 3148: return 255.091084173;    // 255Md
		case 3149: return 256.0938913;      // 256Md
		case 3150: return 257.095542429;    // 257Md
		case 3151: return 258.098431550;    // 258Md
		case 3152: return 259.1005122;      // 259Md
		case 3153: return 260.1036534;      // 260Md
		case 3154: return 261.1058362;      // 261Md
		case 3155: return 262.1091045;      // 262Md
		case 3156: return 248.0865524;      // 248No
		case 3157: return 249.0878030;      // 249No
		case 3158: return 250.0875622;      // 250No
		case 3159: return 251.0889412;      // 251No
		case 3160: return 252.08896710;     // 252No
		case 3161: return 253.090564175;    // 253No
		case 3162: return 254.09095611;     // 254No
		case 3163: return 255.09319116;     // 255No
		case 3164: return 256.094282984;    // 256No
		case 3165: return 257.096887874;    // 257No
		case 3166: return 258.0982111;      // 258No
		case 3167: return 259.1010311;      // 259No
		case 3168: return 260.1026422;      // 260No
		case 3169: return 261.1057022;      // 261No
		case 3170: return 262.1074639;      // 262No
		case 3171: return 263.1107153;      // 263No
		case 3172: return 264.1127370;      // 264No
		case 3173: return 251.0941832;      // 251Lr
		case 3174: return 252.0952626;      // 252Lr
		case 3175: return 253.0950922;      // 253Lr
		case 3176: return 254.0964832;      // 254Lr
		case 3177: return 255.09656219;     // 255Lr
		case 3178: return 256.09849489;     // 256Lr
		case 3179: return 257.09941847;     // 257Lr
		case 3180: return 258.1017611;      // 258Lr
		case 3181: return 259.10290276;     // 259Lr
		case 3182: return 260.1055013;      // 260Lr
		case 3183: return 261.1068822;      // 261Lr
		case 3184: return 262.1096122;      // 262Lr
		case 3185: return 263.1113630;      // 263Lr
		case 3186: return 264.1142047;      // 264Lr
		case 3187: return 265.1161965;      // 265Lr
		case 3188: return 266.1198356;      // 266Lr
		case 3189: return 253.1004444;      // 253Rf
		case 3190: return 254.1000530;      // 254Rf
		case 3191: return 255.1012712;      // 255Rf
		case 3192: return 256.10115219;     // 256Rf
		case 3193: return 257.10291812;     // 257Rf
		case 3194: return 258.10342834;     // 258Rf
		case 3195: return 259.10559678;     // 259Rf
		case 3196: return 260.1064422;      // 260Rf
		case 3197: return 261.10877354;     // 261Rf
		case 3198: return 262.1099224;      // 262Rf
		case 3199: return 263.1124920;      // 263Rf
		case 3200: return 264.1138839;      // 264Rf
		case 3201: return 265.1166839;      // 265Rf
		case 3202: return 266.1181750;      // 266Rf
		case 3203: return 267.1217962;      // 267Rf
		case 3204: return 268.1239777;      // 268Rf
		case 3205: return 255.1070745;      // 255Db
		case 3206: return 256.1078926;      // 256Db
		case 3207: return 257.1075822;      // 257Db
		case 3208: return 258.1092833;      // 258Db
		case 3209: return 259.10949257;     // 259Db
		case 3210: return 260.1113010;      // 260Db
		case 3211: return 261.1119212;      // 261Db
		case 3212: return 262.1140715;      // 262Db
		case 3213: return 263.1149918;      // 263Db
		case 3214: return 264.1174125;      // 264Db
		case 3215: return 265.1186124;      // 265Db
		case 3216: return 266.1210330;      // 266Db
		case 3217: return 267.1224744;      // 267Db
		case 3218: return 268.1256757;      // 268Db
		case 3219: return 269.1279173;      // 269Db
		case 3220: return 270.1313664;      // 270Db
		case 3221: return 258.1129844;      // 258Sg
		case 3222: return 259.1144013;      // 259Sg
		case 3223: return 260.11438422;     // 260Sg
		case 3224: return 261.11594920;     // 261Sg
		case 3225: return 262.11633738;     // 262Sg
		case 3226: return 263.1182910;      // 263Sg
		case 3227: return 264.1189330;      // 264Sg
		case 3228: return 265.1210913;      // 265Sg
		case 3229: return 266.1219826;      // 266Sg
		case 3230: return 267.1243630;      // 267Sg
		case 3231: return 268.1253950;      // 268Sg
		case 3232: return 269.1286339;      // 269Sg
		case 3233: return 270.1304360;      // 270Sg
		case 3234: return 271.1339363;      // 271Sg
		case 3235: return 272.1358983;      // 272Sg
		case 3236: return 273.1395854;      // 273Sg
		case 3237: return 260.1216626;      // 260Bh
		case 3238: return 261.1214522;      // 261Bh
		case 3239: return 262.1229733;      // 262Bh
		case 3240: return 263.1229233;      // 263Bh
		case 3241: return 264.1245919;      // 264Bh
		case 3242: return 265.1249125;      // 265Bh
		case 3243: return 266.1267918;      // 266Bh
		case 3244: return 267.1275028;      // 267Bh
		case 3245: return 268.1296941;      // 268Bh
		case 3246: return 269.1304240;      // 269Bh
		case 3247: return 270.1333631;      // 270Bh
		case 3248: return 271.1352648;      // 271Bh
		case 3249: return 272.1382658;      // 272Bh
		case 3250: return 273.1402480;      // 273Bh
		case 3251: return 274.1435565;      // 274Bh
		case 3252: return 275.1456764;      // 275Bh
		case 3253: return 263.1285214;      // 263Hs
		case 3254: return 264.12835731;     // 264Hs
		case 3255: return 265.12979326;     // 265Hs
		case 3256: return 266.13004642;     // 266Hs
		case 3257: return 267.1316710;      // 267Hs
		case 3258: return 268.1318630;      // 268Hs
		case 3259: return 269.1337513;      // 269Hs
		case 3260: return 270.1342927;      // 270Hs
		case 3261: return 271.1371732;      // 271Hs
		case 3262: return 272.1385055;      // 272Hs
		case 3263: return 273.1416840;      // 273Hs
		case 3264: return 274.1433063;      // 274Hs
		case 3265: return 275.1466763;      // 275Hs
		case 3266: return 276.1484686;      // 276Hs
		case 3267: return 277.1519058;      // 277Hs
		case 3268: return 265.1360048;      // 265Mt
		case 3269: return 266.1373733;      // 266Mt
		case 3270: return 267.1371954;      // 267Mt
		case 3271: return 268.1386525;      // 268Mt
		case 3272: return 269.1388250;      // 269Mt
		case 3273: return 270.1403318;      // 270Mt
		case 3274: return 271.1407435;      // 271Mt
		case 3275: return 272.1434152;      // 272Mt
		case 3276: return 273.1444052;      // 273Mt
		case 3277: return 274.1472438;      // 274Mt
		case 3278: return 275.1488250;      // 275Mt
		case 3279: return 276.1515959;      // 276Mt
		case 3280: return 277.1532782;      // 277Mt
		case 3281: return 278.1563168;      // 278Mt
		case 3282: return 279.1580872;      // 279Mt
		case 3283: return 267.1437715;      // 267Ds
		case 3284: return 268.1434832;      // 268Ds
		case 3285: return 269.14475234;     // 269Ds
		case 3286: return 270.14458452;     // 270Ds
		case 3287: return 271.1459510;      // 271Ds
		case 3288: return 272.1460244;      // 272Ds
		case 3289: return 273.1485614;      // 273Ds
		case 3290: return 274.1494142;      // 274Ds
		case 3291: return 275.1520345;      // 275Ds
		case 3292: return 276.1530359;      // 276Ds
		case 3293: return 277.1559141;      // 277Ds
		case 3294: return 278.1570467;      // 278Ds
		case 3295: return 279.1601064;      // 279Ds
		case 3296: return 280.1613189;      // 280Ds
		case 3297: return 281.1645159;      // 281Ds
		case 3298: return 272.1532725;      // 272Rg
		case 3299: return 273.1531356;      // 273Rg
		case 3300: return 274.1552519;      // 274Rg
		case 3301: return 275.1559456;      // 275Rg
		case 3302: return 276.1583368;      // 276Rg
		case 3303: return 277.1590761;      // 277Rg
		case 3304: return 278.1614938;      // 278Rg
		case 3305: return 279.1627251;      // 279Rg
		case 3306: return 280.1651461;      // 280Rg
		case 3307: return 281.1663689;      // 281Rg
		case 3308: return 282.1691272;      // 282Rg
		case 3309: return 283.1705479;      // 283Rg
		case 3310: return 276.1614164;      // 276Cn
		case 3311: return 277.1636415;      // 277Cn
		case 3312: return 278.1641647;      // 278Cn
		case 3313: return 279.1665450;      // 279Cn
		case 3314: return 280.1671563;      // 280Cn
		case 3315: return 281.1697542;      // 281Cn
		case 3316: return 282.1705070;      // 282Cn
		case 3317: return 283.1732765;      // 283Cn
		case 3318: return 284.1741691;      // 284Cn
		case 3319: return 285.1771260;      // 285Cn
		case 3320: return 278.1705820;      // 278Nh
		case 3321: return 279.1709575;      // 279Nh
		case 3322: return 280.1729375;      // 280Nh
		case 3323: return 281.1734875;      // 281Nh
		case 3324: return 282.1756739;      // 282Nh
		case 3325: return 283.1765752;      // 283Nh
		case 3326: return 284.1787362;      // 284Nh
		case 3327: return 285.1797389;      // 285Nh
		case 3328: return 286.1822172;      // 286Nh
		case 3329: return 287.1833981;      // 287Nh
		case 3330: return 285.1836447;      // 285Fl
		case 3331: return 286.1842371;      // 286Fl
		case 3332: return 287.1867866;      // 287Fl
		case 3333: return 288.1875791;      // 288Fl
		case 3334: return 289.1904260;      // 289Fl
		case 3335: return 287.1907052;      // 287Mc
		case 3336: return 288.1927462;      // 288Mc
		case 3337: return 289.1936389;      // 289Mc
		case 3338: return 290.1959873;      // 290Mc
		case 3339: return 291.1970788;      // 291Uup
		case 3340: return 289.1981657;      // 289Lv
		case 3341: return 290.1986471;      // 290Lv
		case 3342: return 291.2010866;      // 291Lv
		case 3343: return 292.2017491;      // 292Lv
		case 3344: return 293.2044960;      // 293Lv
		case 3345: return 291.2055368;      // 291Ts
		case 3346: return 292.2074675;      // 292Ts
		case 3347: return 293.2082489;      // 293Ts
		case 3348: return 294.2104674;      // 294Uus
		case 3349: return 293.2135678;      // 293Og
		case 3350: return 294.2139271;      // 294Og
		case 3351: return 295.2162469;      // 295Og
		  default: return 0.0;
	}
}

/******************************************************************************

 Function nist_atomic_symbol(): return the atomic symbol for a given isotope.
 See the isotope enumeration in the header file for details.

******************************************************************************/

const char *nist_atomic_symbol(const isotope a)
{
	switch (a)
	{
		case    0: return "1H";
		case    1: return "2H";
		case    2: return "3H";
		case    3: return "4H";
		case    4: return "5H";
		case    5: return "6H";
		case    6: return "7H";
		case    7: return "3He";
		case    8: return "4He";
		case    9: return "5He";
		case   10: return "6He";
		case   11: return "7He";
		case   12: return "8He";
		case   13: return "9He";
		case   14: return "10He";
		case   15: return "3Li";
		case   16: return "4Li";
		case   17: return "5Li";
		case   18: return "6Li";
		case   19: return "7Li";
		case   20: return "8Li";
		case   21: return "9Li";
		case   22: return "10Li";
		case   23: return "11Li";
		case   24: return "12Li";
		case   25: return "13Li";
		case   26: return "5Be";
		case   27: return "6Be";
		case   28: return "7Be";
		case   29: return "8Be";
		case   30: return "9Be";
		case   31: return "10Be";
		case   32: return "11Be";
		case   33: return "12Be";
		case   34: return "13Be";
		case   35: return "14Be";
		case   36: return "15Be";
		case   37: return "16Be";
		case   38: return "6B";
		case   39: return "7B";
		case   40: return "8B";
		case   41: return "9B";
		case   42: return "10B";
		case   43: return "11B";
		case   44: return "12B";
		case   45: return "13B";
		case   46: return "14B";
		case   47: return "15B";
		case   48: return "16B";
		case   49: return "17B";
		case   50: return "18B";
		case   51: return "19B";
		case   52: return "20B";
		case   53: return "21B";
		case   54: return "8C";
		case   55: return "9C";
		case   56: return "10C";
		case   57: return "11C";
		case   58: return "12C";
		case   59: return "13C";
		case   60: return "14C";
		case   61: return "15C";
		case   62: return "16C";
		case   63: return "17C";
		case   64: return "18C";
		case   65: return "19C";
		case   66: return "20C";
		case   67: return "21C";
		case   68: return "22C";
		case   69: return "23C";
		case   70: return "10N";
		case   71: return "11N";
		case   72: return "12N";
		case   73: return "13N";
		case   74: return "14N";
		case   75: return "15N";
		case   76: return "16N";
		case   77: return "17N";
		case   78: return "18N";
		case   79: return "19N";
		case   80: return "20N";
		case   81: return "21N";
		case   82: return "22N";
		case   83: return "23N";
		case   84: return "24N";
		case   85: return "25N";
		case   86: return "12O";
		case   87: return "13O";
		case   88: return "14O";
		case   89: return "15O";
		case   90: return "16O";
		case   91: return "17O";
		case   92: return "18O";
		case   93: return "19O";
		case   94: return "20O";
		case   95: return "21O";
		case   96: return "22O";
		case   97: return "23O";
		case   98: return "24O";
		case   99: return "25O";
		case  100: return "26O";
		case  101: return "27O";
		case  102: return "28O";
		case  103: return "14F";
		case  104: return "15F";
		case  105: return "16F";
		case  106: return "17F";
		case  107: return "18F";
		case  108: return "19F";
		case  109: return "20F";
		case  110: return "21F";
		case  111: return "22F";
		case  112: return "23F";
		case  113: return "24F";
		case  114: return "25F";
		case  115: return "26F";
		case  116: return "27F";
		case  117: return "28F";
		case  118: return "29F";
		case  119: return "30F";
		case  120: return "31F";
		case  121: return "16Ne";
		case  122: return "17Ne";
		case  123: return "18Ne";
		case  124: return "19Ne";
		case  125: return "20Ne";
		case  126: return "21Ne";
		case  127: return "22Ne";
		case  128: return "23Ne";
		case  129: return "24Ne";
		case  130: return "25Ne";
		case  131: return "26Ne";
		case  132: return "27Ne";
		case  133: return "28Ne";
		case  134: return "29Ne";
		case  135: return "30Ne";
		case  136: return "31Ne";
		case  137: return "32Ne";
		case  138: return "33Ne";
		case  139: return "34Ne";
		case  140: return "18Na";
		case  141: return "19Na";
		case  142: return "20Na";
		case  143: return "21Na";
		case  144: return "22Na";
		case  145: return "23Na";
		case  146: return "24Na";
		case  147: return "25Na";
		case  148: return "26Na";
		case  149: return "27Na";
		case  150: return "28Na";
		case  151: return "29Na";
		case  152: return "30Na";
		case  153: return "31Na";
		case  154: return "32Na";
		case  155: return "33Na";
		case  156: return "34Na";
		case  157: return "35Na";
		case  158: return "36Na";
		case  159: return "37Na";
		case  160: return "19Mg";
		case  161: return "20Mg";
		case  162: return "21Mg";
		case  163: return "22Mg";
		case  164: return "23Mg";
		case  165: return "24Mg";
		case  166: return "25Mg";
		case  167: return "26Mg";
		case  168: return "27Mg";
		case  169: return "28Mg";
		case  170: return "29Mg";
		case  171: return "30Mg";
		case  172: return "31Mg";
		case  173: return "32Mg";
		case  174: return "33Mg";
		case  175: return "34Mg";
		case  176: return "35Mg";
		case  177: return "36Mg";
		case  178: return "37Mg";
		case  179: return "38Mg";
		case  180: return "39Mg";
		case  181: return "40Mg";
		case  182: return "21Al";
		case  183: return "22Al";
		case  184: return "23Al";
		case  185: return "24Al";
		case  186: return "25Al";
		case  187: return "26Al";
		case  188: return "27Al";
		case  189: return "28Al";
		case  190: return "29Al";
		case  191: return "30Al";
		case  192: return "31Al";
		case  193: return "32Al";
		case  194: return "33Al";
		case  195: return "34Al";
		case  196: return "35Al";
		case  197: return "36Al";
		case  198: return "37Al";
		case  199: return "38Al";
		case  200: return "39Al";
		case  201: return "40Al";
		case  202: return "41Al";
		case  203: return "42Al";
		case  204: return "43Al";
		case  205: return "22Si";
		case  206: return "23Si";
		case  207: return "24Si";
		case  208: return "25Si";
		case  209: return "26Si";
		case  210: return "27Si";
		case  211: return "28Si";
		case  212: return "29Si";
		case  213: return "30Si";
		case  214: return "31Si";
		case  215: return "32Si";
		case  216: return "33Si";
		case  217: return "34Si";
		case  218: return "35Si";
		case  219: return "36Si";
		case  220: return "37Si";
		case  221: return "38Si";
		case  222: return "39Si";
		case  223: return "40Si";
		case  224: return "41Si";
		case  225: return "42Si";
		case  226: return "43Si";
		case  227: return "44Si";
		case  228: return "45Si";
		case  229: return "24P";
		case  230: return "25P";
		case  231: return "26P";
		case  232: return "27P";
		case  233: return "28P";
		case  234: return "29P";
		case  235: return "30P";
		case  236: return "31P";
		case  237: return "32P";
		case  238: return "33P";
		case  239: return "34P";
		case  240: return "35P";
		case  241: return "36P";
		case  242: return "37P";
		case  243: return "38P";
		case  244: return "39P";
		case  245: return "40P";
		case  246: return "41P";
		case  247: return "42P";
		case  248: return "43P";
		case  249: return "44P";
		case  250: return "45P";
		case  251: return "46P";
		case  252: return "47P";
		case  253: return "26S";
		case  254: return "27S";
		case  255: return "28S";
		case  256: return "29S";
		case  257: return "30S";
		case  258: return "31S";
		case  259: return "32S";
		case  260: return "33S";
		case  261: return "34S";
		case  262: return "35S";
		case  263: return "36S";
		case  264: return "37S";
		case  265: return "38S";
		case  266: return "39S";
		case  267: return "40S";
		case  268: return "41S";
		case  269: return "42S";
		case  270: return "43S";
		case  271: return "44S";
		case  272: return "45S";
		case  273: return "46S";
		case  274: return "47S";
		case  275: return "48S";
		case  276: return "49S";
		case  277: return "28Cl";
		case  278: return "29Cl";
		case  279: return "30Cl";
		case  280: return "31Cl";
		case  281: return "32Cl";
		case  282: return "33Cl";
		case  283: return "34Cl";
		case  284: return "35Cl";
		case  285: return "36Cl";
		case  286: return "37Cl";
		case  287: return "38Cl";
		case  288: return "39Cl";
		case  289: return "40Cl";
		case  290: return "41Cl";
		case  291: return "42Cl";
		case  292: return "43Cl";
		case  293: return "44Cl";
		case  294: return "45Cl";
		case  295: return "46Cl";
		case  296: return "47Cl";
		case  297: return "48Cl";
		case  298: return "49Cl";
		case  299: return "50Cl";
		case  300: return "51Cl";
		case  301: return "30Ar";
		case  302: return "31Ar";
		case  303: return "32Ar";
		case  304: return "33Ar";
		case  305: return "34Ar";
		case  306: return "35Ar";
		case  307: return "36Ar";
		case  308: return "37Ar";
		case  309: return "38Ar";
		case  310: return "39Ar";
		case  311: return "40Ar";
		case  312: return "41Ar";
		case  313: return "42Ar";
		case  314: return "43Ar";
		case  315: return "44Ar";
		case  316: return "45Ar";
		case  317: return "46Ar";
		case  318: return "47Ar";
		case  319: return "48Ar";
		case  320: return "49Ar";
		case  321: return "50Ar";
		case  322: return "51Ar";
		case  323: return "52Ar";
		case  324: return "53Ar";
		case  325: return "32K";
		case  326: return "33K";
		case  327: return "34K";
		case  328: return "35K";
		case  329: return "36K";
		case  330: return "37K";
		case  331: return "38K";
		case  332: return "39K";
		case  333: return "40K";
		case  334: return "41K";
		case  335: return "42K";
		case  336: return "43K";
		case  337: return "44K";
		case  338: return "45K";
		case  339: return "46K";
		case  340: return "47K";
		case  341: return "48K";
		case  342: return "49K";
		case  343: return "50K";
		case  344: return "51K";
		case  345: return "52K";
		case  346: return "53K";
		case  347: return "54K";
		case  348: return "55K";
		case  349: return "56K";
		case  350: return "34Ca";
		case  351: return "35Ca";
		case  352: return "36Ca";
		case  353: return "37Ca";
		case  354: return "38Ca";
		case  355: return "39Ca";
		case  356: return "40Ca";
		case  357: return "41Ca";
		case  358: return "42Ca";
		case  359: return "43Ca";
		case  360: return "44Ca";
		case  361: return "45Ca";
		case  362: return "46Ca";
		case  363: return "47Ca";
		case  364: return "48Ca";
		case  365: return "49Ca";
		case  366: return "50Ca";
		case  367: return "51Ca";
		case  368: return "52Ca";
		case  369: return "53Ca";
		case  370: return "54Ca";
		case  371: return "55Ca";
		case  372: return "56Ca";
		case  373: return "57Ca";
		case  374: return "58Ca";
		case  375: return "36Sc";
		case  376: return "37Sc";
		case  377: return "38Sc";
		case  378: return "39Sc";
		case  379: return "40Sc";
		case  380: return "41Sc";
		case  381: return "42Sc";
		case  382: return "43Sc";
		case  383: return "44Sc";
		case  384: return "45Sc";
		case  385: return "46Sc";
		case  386: return "47Sc";
		case  387: return "48Sc";
		case  388: return "49Sc";
		case  389: return "50Sc";
		case  390: return "51Sc";
		case  391: return "52Sc";
		case  392: return "53Sc";
		case  393: return "54Sc";
		case  394: return "55Sc";
		case  395: return "56Sc";
		case  396: return "57Sc";
		case  397: return "58Sc";
		case  398: return "59Sc";
		case  399: return "60Sc";
		case  400: return "61Sc";
		case  401: return "38Ti";
		case  402: return "39Ti";
		case  403: return "40Ti";
		case  404: return "41Ti";
		case  405: return "42Ti";
		case  406: return "43Ti";
		case  407: return "44Ti";
		case  408: return "45Ti";
		case  409: return "46Ti";
		case  410: return "47Ti";
		case  411: return "48Ti";
		case  412: return "49Ti";
		case  413: return "50Ti";
		case  414: return "51Ti";
		case  415: return "52Ti";
		case  416: return "53Ti";
		case  417: return "54Ti";
		case  418: return "55Ti";
		case  419: return "56Ti";
		case  420: return "57Ti";
		case  421: return "58Ti";
		case  422: return "59Ti";
		case  423: return "60Ti";
		case  424: return "61Ti";
		case  425: return "62Ti";
		case  426: return "63Ti";
		case  427: return "40V";
		case  428: return "41V";
		case  429: return "42V";
		case  430: return "43V";
		case  431: return "44V";
		case  432: return "45V";
		case  433: return "46V";
		case  434: return "47V";
		case  435: return "48V";
		case  436: return "49V";
		case  437: return "50V";
		case  438: return "51V";
		case  439: return "52V";
		case  440: return "53V";
		case  441: return "54V";
		case  442: return "55V";
		case  443: return "56V";
		case  444: return "57V";
		case  445: return "58V";
		case  446: return "59V";
		case  447: return "60V";
		case  448: return "61V";
		case  449: return "62V";
		case  450: return "63V";
		case  451: return "64V";
		case  452: return "65V";
		case  453: return "66V";
		case  454: return "42Cr";
		case  455: return "43Cr";
		case  456: return "44Cr";
		case  457: return "45Cr";
		case  458: return "46Cr";
		case  459: return "47Cr";
		case  460: return "48Cr";
		case  461: return "49Cr";
		case  462: return "50Cr";
		case  463: return "51Cr";
		case  464: return "52Cr";
		case  465: return "53Cr";
		case  466: return "54Cr";
		case  467: return "55Cr";
		case  468: return "56Cr";
		case  469: return "57Cr";
		case  470: return "58Cr";
		case  471: return "59Cr";
		case  472: return "60Cr";
		case  473: return "61Cr";
		case  474: return "62Cr";
		case  475: return "63Cr";
		case  476: return "64Cr";
		case  477: return "65Cr";
		case  478: return "66Cr";
		case  479: return "67Cr";
		case  480: return "68Cr";
		case  481: return "44Mn";
		case  482: return "45Mn";
		case  483: return "46Mn";
		case  484: return "47Mn";
		case  485: return "48Mn";
		case  486: return "49Mn";
		case  487: return "50Mn";
		case  488: return "51Mn";
		case  489: return "52Mn";
		case  490: return "53Mn";
		case  491: return "54Mn";
		case  492: return "55Mn";
		case  493: return "56Mn";
		case  494: return "57Mn";
		case  495: return "58Mn";
		case  496: return "59Mn";
		case  497: return "60Mn";
		case  498: return "61Mn";
		case  499: return "62Mn";
		case  500: return "63Mn";
		case  501: return "64Mn";
		case  502: return "65Mn";
		case  503: return "66Mn";
		case  504: return "67Mn";
		case  505: return "68Mn";
		case  506: return "69Mn";
		case  507: return "70Mn";
		case  508: return "71Mn";
		case  509: return "45Fe";
		case  510: return "46Fe";
		case  511: return "47Fe";
		case  512: return "48Fe";
		case  513: return "49Fe";
		case  514: return "50Fe";
		case  515: return "51Fe";
		case  516: return "52Fe";
		case  517: return "53Fe";
		case  518: return "54Fe";
		case  519: return "55Fe";
		case  520: return "56Fe";
		case  521: return "57Fe";
		case  522: return "58Fe";
		case  523: return "59Fe";
		case  524: return "60Fe";
		case  525: return "61Fe";
		case  526: return "62Fe";
		case  527: return "63Fe";
		case  528: return "64Fe";
		case  529: return "65Fe";
		case  530: return "66Fe";
		case  531: return "67Fe";
		case  532: return "68Fe";
		case  533: return "69Fe";
		case  534: return "70Fe";
		case  535: return "71Fe";
		case  536: return "72Fe";
		case  537: return "73Fe";
		case  538: return "74Fe";
		case  539: return "47Co";
		case  540: return "48Co";
		case  541: return "49Co";
		case  542: return "50Co";
		case  543: return "51Co";
		case  544: return "52Co";
		case  545: return "53Co";
		case  546: return "54Co";
		case  547: return "55Co";
		case  548: return "56Co";
		case  549: return "57Co";
		case  550: return "58Co";
		case  551: return "59Co";
		case  552: return "60Co";
		case  553: return "61Co";
		case  554: return "62Co";
		case  555: return "63Co";
		case  556: return "64Co";
		case  557: return "65Co";
		case  558: return "66Co";
		case  559: return "67Co";
		case  560: return "68Co";
		case  561: return "69Co";
		case  562: return "70Co";
		case  563: return "71Co";
		case  564: return "72Co";
		case  565: return "73Co";
		case  566: return "74Co";
		case  567: return "75Co";
		case  568: return "76Co";
		case  569: return "48Ni";
		case  570: return "49Ni";
		case  571: return "50Ni";
		case  572: return "51Ni";
		case  573: return "52Ni";
		case  574: return "53Ni";
		case  575: return "54Ni";
		case  576: return "55Ni";
		case  577: return "56Ni";
		case  578: return "57Ni";
		case  579: return "58Ni";
		case  580: return "59Ni";
		case  581: return "60Ni";
		case  582: return "61Ni";
		case  583: return "62Ni";
		case  584: return "63Ni";
		case  585: return "64Ni";
		case  586: return "65Ni";
		case  587: return "66Ni";
		case  588: return "67Ni";
		case  589: return "68Ni";
		case  590: return "69Ni";
		case  591: return "70Ni";
		case  592: return "71Ni";
		case  593: return "72Ni";
		case  594: return "73Ni";
		case  595: return "74Ni";
		case  596: return "75Ni";
		case  597: return "76Ni";
		case  598: return "77Ni";
		case  599: return "78Ni";
		case  600: return "79Ni";
		case  601: return "52Cu";
		case  602: return "53Cu";
		case  603: return "54Cu";
		case  604: return "55Cu";
		case  605: return "56Cu";
		case  606: return "57Cu";
		case  607: return "58Cu";
		case  608: return "59Cu";
		case  609: return "60Cu";
		case  610: return "61Cu";
		case  611: return "62Cu";
		case  612: return "63Cu";
		case  613: return "64Cu";
		case  614: return "65Cu";
		case  615: return "66Cu";
		case  616: return "67Cu";
		case  617: return "68Cu";
		case  618: return "69Cu";
		case  619: return "70Cu";
		case  620: return "71Cu";
		case  621: return "72Cu";
		case  622: return "73Cu";
		case  623: return "74Cu";
		case  624: return "75Cu";
		case  625: return "76Cu";
		case  626: return "77Cu";
		case  627: return "78Cu";
		case  628: return "79Cu";
		case  629: return "80Cu";
		case  630: return "81Cu";
		case  631: return "82Cu";
		case  632: return "54Zn";
		case  633: return "55Zn";
		case  634: return "56Zn";
		case  635: return "57Zn";
		case  636: return "58Zn";
		case  637: return "59Zn";
		case  638: return "60Zn";
		case  639: return "61Zn";
		case  640: return "62Zn";
		case  641: return "63Zn";
		case  642: return "64Zn";
		case  643: return "65Zn";
		case  644: return "66Zn";
		case  645: return "67Zn";
		case  646: return "68Zn";
		case  647: return "69Zn";
		case  648: return "70Zn";
		case  649: return "71Zn";
		case  650: return "72Zn";
		case  651: return "73Zn";
		case  652: return "74Zn";
		case  653: return "75Zn";
		case  654: return "76Zn";
		case  655: return "77Zn";
		case  656: return "78Zn";
		case  657: return "79Zn";
		case  658: return "80Zn";
		case  659: return "81Zn";
		case  660: return "82Zn";
		case  661: return "83Zn";
		case  662: return "84Zn";
		case  663: return "85Zn";
		case  664: return "56Ga";
		case  665: return "57Ga";
		case  666: return "58Ga";
		case  667: return "59Ga";
		case  668: return "60Ga";
		case  669: return "61Ga";
		case  670: return "62Ga";
		case  671: return "63Ga";
		case  672: return "64Ga";
		case  673: return "65Ga";
		case  674: return "66Ga";
		case  675: return "67Ga";
		case  676: return "68Ga";
		case  677: return "69Ga";
		case  678: return "70Ga";
		case  679: return "71Ga";
		case  680: return "72Ga";
		case  681: return "73Ga";
		case  682: return "74Ga";
		case  683: return "75Ga";
		case  684: return "76Ga";
		case  685: return "77Ga";
		case  686: return "78Ga";
		case  687: return "79Ga";
		case  688: return "80Ga";
		case  689: return "81Ga";
		case  690: return "82Ga";
		case  691: return "83Ga";
		case  692: return "84Ga";
		case  693: return "85Ga";
		case  694: return "86Ga";
		case  695: return "87Ga";
		case  696: return "58Ge";
		case  697: return "59Ge";
		case  698: return "60Ge";
		case  699: return "61Ge";
		case  700: return "62Ge";
		case  701: return "63Ge";
		case  702: return "64Ge";
		case  703: return "65Ge";
		case  704: return "66Ge";
		case  705: return "67Ge";
		case  706: return "68Ge";
		case  707: return "69Ge";
		case  708: return "70Ge";
		case  709: return "71Ge";
		case  710: return "72Ge";
		case  711: return "73Ge";
		case  712: return "74Ge";
		case  713: return "75Ge";
		case  714: return "76Ge";
		case  715: return "77Ge";
		case  716: return "78Ge";
		case  717: return "79Ge";
		case  718: return "80Ge";
		case  719: return "81Ge";
		case  720: return "82Ge";
		case  721: return "83Ge";
		case  722: return "84Ge";
		case  723: return "85Ge";
		case  724: return "86Ge";
		case  725: return "87Ge";
		case  726: return "88Ge";
		case  727: return "89Ge";
		case  728: return "90Ge";
		case  729: return "60As";
		case  730: return "61As";
		case  731: return "62As";
		case  732: return "63As";
		case  733: return "64As";
		case  734: return "65As";
		case  735: return "66As";
		case  736: return "67As";
		case  737: return "68As";
		case  738: return "69As";
		case  739: return "70As";
		case  740: return "71As";
		case  741: return "72As";
		case  742: return "73As";
		case  743: return "74As";
		case  744: return "75As";
		case  745: return "76As";
		case  746: return "77As";
		case  747: return "78As";
		case  748: return "79As";
		case  749: return "80As";
		case  750: return "81As";
		case  751: return "82As";
		case  752: return "83As";
		case  753: return "84As";
		case  754: return "85As";
		case  755: return "86As";
		case  756: return "87As";
		case  757: return "88As";
		case  758: return "89As";
		case  759: return "90As";
		case  760: return "91As";
		case  761: return "92As";
		case  762: return "64Se";
		case  763: return "65Se";
		case  764: return "66Se";
		case  765: return "67Se";
		case  766: return "68Se";
		case  767: return "69Se";
		case  768: return "70Se";
		case  769: return "71Se";
		case  770: return "72Se";
		case  771: return "73Se";
		case  772: return "74Se";
		case  773: return "75Se";
		case  774: return "76Se";
		case  775: return "77Se";
		case  776: return "78Se";
		case  777: return "79Se";
		case  778: return "80Se";
		case  779: return "81Se";
		case  780: return "82Se";
		case  781: return "83Se";
		case  782: return "84Se";
		case  783: return "85Se";
		case  784: return "86Se";
		case  785: return "87Se";
		case  786: return "88Se";
		case  787: return "89Se";
		case  788: return "90Se";
		case  789: return "91Se";
		case  790: return "92Se";
		case  791: return "93Se";
		case  792: return "94Se";
		case  793: return "95Se";
		case  794: return "67Br";
		case  795: return "68Br";
		case  796: return "69Br";
		case  797: return "70Br";
		case  798: return "71Br";
		case  799: return "72Br";
		case  800: return "73Br";
		case  801: return "74Br";
		case  802: return "75Br";
		case  803: return "76Br";
		case  804: return "77Br";
		case  805: return "78Br";
		case  806: return "79Br";
		case  807: return "80Br";
		case  808: return "81Br";
		case  809: return "82Br";
		case  810: return "83Br";
		case  811: return "84Br";
		case  812: return "85Br";
		case  813: return "86Br";
		case  814: return "87Br";
		case  815: return "88Br";
		case  816: return "89Br";
		case  817: return "90Br";
		case  818: return "91Br";
		case  819: return "92Br";
		case  820: return "93Br";
		case  821: return "94Br";
		case  822: return "95Br";
		case  823: return "96Br";
		case  824: return "97Br";
		case  825: return "98Br";
		case  826: return "69Kr";
		case  827: return "70Kr";
		case  828: return "71Kr";
		case  829: return "72Kr";
		case  830: return "73Kr";
		case  831: return "74Kr";
		case  832: return "75Kr";
		case  833: return "76Kr";
		case  834: return "77Kr";
		case  835: return "78Kr";
		case  836: return "79Kr";
		case  837: return "80Kr";
		case  838: return "81Kr";
		case  839: return "82Kr";
		case  840: return "83Kr";
		case  841: return "84Kr";
		case  842: return "85Kr";
		case  843: return "86Kr";
		case  844: return "87Kr";
		case  845: return "88Kr";
		case  846: return "89Kr";
		case  847: return "90Kr";
		case  848: return "91Kr";
		case  849: return "92Kr";
		case  850: return "93Kr";
		case  851: return "94Kr";
		case  852: return "95Kr";
		case  853: return "96Kr";
		case  854: return "97Kr";
		case  855: return "98Kr";
		case  856: return "99Kr";
		case  857: return "100Kr";
		case  858: return "101Kr";
		case  859: return "71Rb";
		case  860: return "72Rb";
		case  861: return "73Rb";
		case  862: return "74Rb";
		case  863: return "75Rb";
		case  864: return "76Rb";
		case  865: return "77Rb";
		case  866: return "78Rb";
		case  867: return "79Rb";
		case  868: return "80Rb";
		case  869: return "81Rb";
		case  870: return "82Rb";
		case  871: return "83Rb";
		case  872: return "84Rb";
		case  873: return "85Rb";
		case  874: return "86Rb";
		case  875: return "87Rb";
		case  876: return "88Rb";
		case  877: return "89Rb";
		case  878: return "90Rb";
		case  879: return "91Rb";
		case  880: return "92Rb";
		case  881: return "93Rb";
		case  882: return "94Rb";
		case  883: return "95Rb";
		case  884: return "96Rb";
		case  885: return "97Rb";
		case  886: return "98Rb";
		case  887: return "99Rb";
		case  888: return "100Rb";
		case  889: return "101Rb";
		case  890: return "102Rb";
		case  891: return "103Rb";
		case  892: return "73Sr";
		case  893: return "74Sr";
		case  894: return "75Sr";
		case  895: return "76Sr";
		case  896: return "77Sr";
		case  897: return "78Sr";
		case  898: return "79Sr";
		case  899: return "80Sr";
		case  900: return "81Sr";
		case  901: return "82Sr";
		case  902: return "83Sr";
		case  903: return "84Sr";
		case  904: return "85Sr";
		case  905: return "86Sr";
		case  906: return "87Sr";
		case  907: return "88Sr";
		case  908: return "89Sr";
		case  909: return "90Sr";
		case  910: return "91Sr";
		case  911: return "92Sr";
		case  912: return "93Sr";
		case  913: return "94Sr";
		case  914: return "95Sr";
		case  915: return "96Sr";
		case  916: return "97Sr";
		case  917: return "98Sr";
		case  918: return "99Sr";
		case  919: return "100Sr";
		case  920: return "101Sr";
		case  921: return "102Sr";
		case  922: return "103Sr";
		case  923: return "104Sr";
		case  924: return "105Sr";
		case  925: return "106Sr";
		case  926: return "107Sr";
		case  927: return "76Y";
		case  928: return "77Y";
		case  929: return "78Y";
		case  930: return "79Y";
		case  931: return "80Y";
		case  932: return "81Y";
		case  933: return "82Y";
		case  934: return "83Y";
		case  935: return "84Y";
		case  936: return "85Y";
		case  937: return "86Y";
		case  938: return "87Y";
		case  939: return "88Y";
		case  940: return "89Y";
		case  941: return "90Y";
		case  942: return "91Y";
		case  943: return "92Y";
		case  944: return "93Y";
		case  945: return "94Y";
		case  946: return "95Y";
		case  947: return "96Y";
		case  948: return "97Y";
		case  949: return "98Y";
		case  950: return "99Y";
		case  951: return "100Y";
		case  952: return "101Y";
		case  953: return "102Y";
		case  954: return "103Y";
		case  955: return "104Y";
		case  956: return "105Y";
		case  957: return "106Y";
		case  958: return "107Y";
		case  959: return "108Y";
		case  960: return "109Y";
		case  961: return "78Zr";
		case  962: return "79Zr";
		case  963: return "80Zr";
		case  964: return "81Zr";
		case  965: return "82Zr";
		case  966: return "83Zr";
		case  967: return "84Zr";
		case  968: return "85Zr";
		case  969: return "86Zr";
		case  970: return "87Zr";
		case  971: return "88Zr";
		case  972: return "89Zr";
		case  973: return "90Zr";
		case  974: return "91Zr";
		case  975: return "92Zr";
		case  976: return "93Zr";
		case  977: return "94Zr";
		case  978: return "95Zr";
		case  979: return "96Zr";
		case  980: return "97Zr";
		case  981: return "98Zr";
		case  982: return "99Zr";
		case  983: return "100Zr";
		case  984: return "101Zr";
		case  985: return "102Zr";
		case  986: return "103Zr";
		case  987: return "104Zr";
		case  988: return "105Zr";
		case  989: return "106Zr";
		case  990: return "107Zr";
		case  991: return "108Zr";
		case  992: return "109Zr";
		case  993: return "110Zr";
		case  994: return "111Zr";
		case  995: return "112Zr";
		case  996: return "81Nb";
		case  997: return "82Nb";
		case  998: return "83Nb";
		case  999: return "84Nb";
		case 1000: return "85Nb";
		case 1001: return "86Nb";
		case 1002: return "87Nb";
		case 1003: return "88Nb";
		case 1004: return "89Nb";
		case 1005: return "90Nb";
		case 1006: return "91Nb";
		case 1007: return "92Nb";
		case 1008: return "93Nb";
		case 1009: return "94Nb";
		case 1010: return "95Nb";
		case 1011: return "96Nb";
		case 1012: return "97Nb";
		case 1013: return "98Nb";
		case 1014: return "99Nb";
		case 1015: return "100Nb";
		case 1016: return "101Nb";
		case 1017: return "102Nb";
		case 1018: return "103Nb";
		case 1019: return "104Nb";
		case 1020: return "105Nb";
		case 1021: return "106Nb";
		case 1022: return "107Nb";
		case 1023: return "108Nb";
		case 1024: return "109Nb";
		case 1025: return "110Nb";
		case 1026: return "111Nb";
		case 1027: return "112Nb";
		case 1028: return "113Nb";
		case 1029: return "114Nb";
		case 1030: return "115Nb";
		case 1031: return "83Mo";
		case 1032: return "84Mo";
		case 1033: return "85Mo";
		case 1034: return "86Mo";
		case 1035: return "87Mo";
		case 1036: return "88Mo";
		case 1037: return "89Mo";
		case 1038: return "90Mo";
		case 1039: return "91Mo";
		case 1040: return "92Mo";
		case 1041: return "93Mo";
		case 1042: return "94Mo";
		case 1043: return "95Mo";
		case 1044: return "96Mo";
		case 1045: return "97Mo";
		case 1046: return "98Mo";
		case 1047: return "99Mo";
		case 1048: return "100Mo";
		case 1049: return "101Mo";
		case 1050: return "102Mo";
		case 1051: return "103Mo";
		case 1052: return "104Mo";
		case 1053: return "105Mo";
		case 1054: return "106Mo";
		case 1055: return "107Mo";
		case 1056: return "108Mo";
		case 1057: return "109Mo";
		case 1058: return "110Mo";
		case 1059: return "111Mo";
		case 1060: return "112Mo";
		case 1061: return "113Mo";
		case 1062: return "114Mo";
		case 1063: return "115Mo";
		case 1064: return "116Mo";
		case 1065: return "117Mo";
		case 1066: return "85Tc";
		case 1067: return "86Tc";
		case 1068: return "87Tc";
		case 1069: return "88Tc";
		case 1070: return "89Tc";
		case 1071: return "90Tc";
		case 1072: return "91Tc";
		case 1073: return "92Tc";
		case 1074: return "93Tc";
		case 1075: return "94Tc";
		case 1076: return "95Tc";
		case 1077: return "96Tc";
		case 1078: return "97Tc";
		case 1079: return "98Tc";
		case 1080: return "99Tc";
		case 1081: return "100Tc";
		case 1082: return "101Tc";
		case 1083: return "102Tc";
		case 1084: return "103Tc";
		case 1085: return "104Tc";
		case 1086: return "105Tc";
		case 1087: return "106Tc";
		case 1088: return "107Tc";
		case 1089: return "108Tc";
		case 1090: return "109Tc";
		case 1091: return "110Tc";
		case 1092: return "111Tc";
		case 1093: return "112Tc";
		case 1094: return "113Tc";
		case 1095: return "114Tc";
		case 1096: return "115Tc";
		case 1097: return "116Tc";
		case 1098: return "117Tc";
		case 1099: return "118Tc";
		case 1100: return "119Tc";
		case 1101: return "120Tc";
		case 1102: return "87Ru";
		case 1103: return "88Ru";
		case 1104: return "89Ru";
		case 1105: return "90Ru";
		case 1106: return "91Ru";
		case 1107: return "92Ru";
		case 1108: return "93Ru";
		case 1109: return "94Ru";
		case 1110: return "95Ru";
		case 1111: return "96Ru";
		case 1112: return "97Ru";
		case 1113: return "98Ru";
		case 1114: return "99Ru";
		case 1115: return "100Ru";
		case 1116: return "101Ru";
		case 1117: return "102Ru";
		case 1118: return "103Ru";
		case 1119: return "104Ru";
		case 1120: return "105Ru";
		case 1121: return "106Ru";
		case 1122: return "107Ru";
		case 1123: return "108Ru";
		case 1124: return "109Ru";
		case 1125: return "110Ru";
		case 1126: return "111Ru";
		case 1127: return "112Ru";
		case 1128: return "113Ru";
		case 1129: return "114Ru";
		case 1130: return "115Ru";
		case 1131: return "116Ru";
		case 1132: return "117Ru";
		case 1133: return "118Ru";
		case 1134: return "119Ru";
		case 1135: return "120Ru";
		case 1136: return "121Ru";
		case 1137: return "122Ru";
		case 1138: return "123Ru";
		case 1139: return "124Ru";
		case 1140: return "89Rh";
		case 1141: return "90Rh";
		case 1142: return "91Rh";
		case 1143: return "92Rh";
		case 1144: return "93Rh";
		case 1145: return "94Rh";
		case 1146: return "95Rh";
		case 1147: return "96Rh";
		case 1148: return "97Rh";
		case 1149: return "98Rh";
		case 1150: return "99Rh";
		case 1151: return "100Rh";
		case 1152: return "101Rh";
		case 1153: return "102Rh";
		case 1154: return "103Rh";
		case 1155: return "104Rh";
		case 1156: return "105Rh";
		case 1157: return "106Rh";
		case 1158: return "107Rh";
		case 1159: return "108Rh";
		case 1160: return "109Rh";
		case 1161: return "110Rh";
		case 1162: return "111Rh";
		case 1163: return "112Rh";
		case 1164: return "113Rh";
		case 1165: return "114Rh";
		case 1166: return "115Rh";
		case 1167: return "116Rh";
		case 1168: return "117Rh";
		case 1169: return "118Rh";
		case 1170: return "119Rh";
		case 1171: return "120Rh";
		case 1172: return "121Rh";
		case 1173: return "122Rh";
		case 1174: return "123Rh";
		case 1175: return "124Rh";
		case 1176: return "125Rh";
		case 1177: return "126Rh";
		case 1178: return "91Pd";
		case 1179: return "92Pd";
		case 1180: return "93Pd";
		case 1181: return "94Pd";
		case 1182: return "95Pd";
		case 1183: return "96Pd";
		case 1184: return "97Pd";
		case 1185: return "98Pd";
		case 1186: return "99Pd";
		case 1187: return "100Pd";
		case 1188: return "101Pd";
		case 1189: return "102Pd";
		case 1190: return "103Pd";
		case 1191: return "104Pd";
		case 1192: return "105Pd";
		case 1193: return "106Pd";
		case 1194: return "107Pd";
		case 1195: return "108Pd";
		case 1196: return "109Pd";
		case 1197: return "110Pd";
		case 1198: return "111Pd";
		case 1199: return "112Pd";
		case 1200: return "113Pd";
		case 1201: return "114Pd";
		case 1202: return "115Pd";
		case 1203: return "116Pd";
		case 1204: return "117Pd";
		case 1205: return "118Pd";
		case 1206: return "119Pd";
		case 1207: return "120Pd";
		case 1208: return "121Pd";
		case 1209: return "122Pd";
		case 1210: return "123Pd";
		case 1211: return "124Pd";
		case 1212: return "125Pd";
		case 1213: return "126Pd";
		case 1214: return "127Pd";
		case 1215: return "128Pd";
		case 1216: return "93Ag";
		case 1217: return "94Ag";
		case 1218: return "95Ag";
		case 1219: return "96Ag";
		case 1220: return "97Ag";
		case 1221: return "98Ag";
		case 1222: return "99Ag";
		case 1223: return "100Ag";
		case 1224: return "101Ag";
		case 1225: return "102Ag";
		case 1226: return "103Ag";
		case 1227: return "104Ag";
		case 1228: return "105Ag";
		case 1229: return "106Ag";
		case 1230: return "107Ag";
		case 1231: return "108Ag";
		case 1232: return "109Ag";
		case 1233: return "110Ag";
		case 1234: return "111Ag";
		case 1235: return "112Ag";
		case 1236: return "113Ag";
		case 1237: return "114Ag";
		case 1238: return "115Ag";
		case 1239: return "116Ag";
		case 1240: return "117Ag";
		case 1241: return "118Ag";
		case 1242: return "119Ag";
		case 1243: return "120Ag";
		case 1244: return "121Ag";
		case 1245: return "122Ag";
		case 1246: return "123Ag";
		case 1247: return "124Ag";
		case 1248: return "125Ag";
		case 1249: return "126Ag";
		case 1250: return "127Ag";
		case 1251: return "128Ag";
		case 1252: return "129Ag";
		case 1253: return "130Ag";
		case 1254: return "95Cd";
		case 1255: return "96Cd";
		case 1256: return "97Cd";
		case 1257: return "98Cd";
		case 1258: return "99Cd";
		case 1259: return "100Cd";
		case 1260: return "101Cd";
		case 1261: return "102Cd";
		case 1262: return "103Cd";
		case 1263: return "104Cd";
		case 1264: return "105Cd";
		case 1265: return "106Cd";
		case 1266: return "107Cd";
		case 1267: return "108Cd";
		case 1268: return "109Cd";
		case 1269: return "110Cd";
		case 1270: return "111Cd";
		case 1271: return "112Cd";
		case 1272: return "113Cd";
		case 1273: return "114Cd";
		case 1274: return "115Cd";
		case 1275: return "116Cd";
		case 1276: return "117Cd";
		case 1277: return "118Cd";
		case 1278: return "119Cd";
		case 1279: return "120Cd";
		case 1280: return "121Cd";
		case 1281: return "122Cd";
		case 1282: return "123Cd";
		case 1283: return "124Cd";
		case 1284: return "125Cd";
		case 1285: return "126Cd";
		case 1286: return "127Cd";
		case 1287: return "128Cd";
		case 1288: return "129Cd";
		case 1289: return "130Cd";
		case 1290: return "131Cd";
		case 1291: return "132Cd";
		case 1292: return "133Cd";
		case 1293: return "97In";
		case 1294: return "98In";
		case 1295: return "99In";
		case 1296: return "100In";
		case 1297: return "101In";
		case 1298: return "102In";
		case 1299: return "103In";
		case 1300: return "104In";
		case 1301: return "105In";
		case 1302: return "106In";
		case 1303: return "107In";
		case 1304: return "108In";
		case 1305: return "109In";
		case 1306: return "110In";
		case 1307: return "111In";
		case 1308: return "112In";
		case 1309: return "113In";
		case 1310: return "114In";
		case 1311: return "115In";
		case 1312: return "116In";
		case 1313: return "117In";
		case 1314: return "118In";
		case 1315: return "119In";
		case 1316: return "120In";
		case 1317: return "121In";
		case 1318: return "122In";
		case 1319: return "123In";
		case 1320: return "124In";
		case 1321: return "125In";
		case 1322: return "126In";
		case 1323: return "127In";
		case 1324: return "128In";
		case 1325: return "129In";
		case 1326: return "130In";
		case 1327: return "131In";
		case 1328: return "132In";
		case 1329: return "133In";
		case 1330: return "134In";
		case 1331: return "135In";
		case 1332: return "99Sn";
		case 1333: return "100Sn";
		case 1334: return "101Sn";
		case 1335: return "102Sn";
		case 1336: return "103Sn";
		case 1337: return "104Sn";
		case 1338: return "105Sn";
		case 1339: return "106Sn";
		case 1340: return "107Sn";
		case 1341: return "108Sn";
		case 1342: return "109Sn";
		case 1343: return "110Sn";
		case 1344: return "111Sn";
		case 1345: return "112Sn";
		case 1346: return "113Sn";
		case 1347: return "114Sn";
		case 1348: return "115Sn";
		case 1349: return "116Sn";
		case 1350: return "117Sn";
		case 1351: return "118Sn";
		case 1352: return "119Sn";
		case 1353: return "120Sn";
		case 1354: return "121Sn";
		case 1355: return "122Sn";
		case 1356: return "123Sn";
		case 1357: return "124Sn";
		case 1358: return "125Sn";
		case 1359: return "126Sn";
		case 1360: return "127Sn";
		case 1361: return "128Sn";
		case 1362: return "129Sn";
		case 1363: return "130Sn";
		case 1364: return "131Sn";
		case 1365: return "132Sn";
		case 1366: return "133Sn";
		case 1367: return "134Sn";
		case 1368: return "135Sn";
		case 1369: return "136Sn";
		case 1370: return "137Sn";
		case 1371: return "138Sn";
		case 1372: return "103Sb";
		case 1373: return "104Sb";
		case 1374: return "105Sb";
		case 1375: return "106Sb";
		case 1376: return "107Sb";
		case 1377: return "108Sb";
		case 1378: return "109Sb";
		case 1379: return "110Sb";
		case 1380: return "111Sb";
		case 1381: return "112Sb";
		case 1382: return "113Sb";
		case 1383: return "114Sb";
		case 1384: return "115Sb";
		case 1385: return "116Sb";
		case 1386: return "117Sb";
		case 1387: return "118Sb";
		case 1388: return "119Sb";
		case 1389: return "120Sb";
		case 1390: return "121Sb";
		case 1391: return "122Sb";
		case 1392: return "123Sb";
		case 1393: return "124Sb";
		case 1394: return "125Sb";
		case 1395: return "126Sb";
		case 1396: return "127Sb";
		case 1397: return "128Sb";
		case 1398: return "129Sb";
		case 1399: return "130Sb";
		case 1400: return "131Sb";
		case 1401: return "132Sb";
		case 1402: return "133Sb";
		case 1403: return "134Sb";
		case 1404: return "135Sb";
		case 1405: return "136Sb";
		case 1406: return "137Sb";
		case 1407: return "138Sb";
		case 1408: return "139Sb";
		case 1409: return "140Sb";
		case 1410: return "105Te";
		case 1411: return "106Te";
		case 1412: return "107Te";
		case 1413: return "108Te";
		case 1414: return "109Te";
		case 1415: return "110Te";
		case 1416: return "111Te";
		case 1417: return "112Te";
		case 1418: return "113Te";
		case 1419: return "114Te";
		case 1420: return "115Te";
		case 1421: return "116Te";
		case 1422: return "117Te";
		case 1423: return "118Te";
		case 1424: return "119Te";
		case 1425: return "120Te";
		case 1426: return "121Te";
		case 1427: return "122Te";
		case 1428: return "123Te";
		case 1429: return "124Te";
		case 1430: return "125Te";
		case 1431: return "126Te";
		case 1432: return "127Te";
		case 1433: return "128Te";
		case 1434: return "129Te";
		case 1435: return "130Te";
		case 1436: return "131Te";
		case 1437: return "132Te";
		case 1438: return "133Te";
		case 1439: return "134Te";
		case 1440: return "135Te";
		case 1441: return "136Te";
		case 1442: return "137Te";
		case 1443: return "138Te";
		case 1444: return "139Te";
		case 1445: return "140Te";
		case 1446: return "141Te";
		case 1447: return "142Te";
		case 1448: return "143Te";
		case 1449: return "107I";
		case 1450: return "108I";
		case 1451: return "109I";
		case 1452: return "110I";
		case 1453: return "111I";
		case 1454: return "112I";
		case 1455: return "113I";
		case 1456: return "114I";
		case 1457: return "115I";
		case 1458: return "116I";
		case 1459: return "117I";
		case 1460: return "118I";
		case 1461: return "119I";
		case 1462: return "120I";
		case 1463: return "121I";
		case 1464: return "122I";
		case 1465: return "123I";
		case 1466: return "124I";
		case 1467: return "125I";
		case 1468: return "126I";
		case 1469: return "127I";
		case 1470: return "128I";
		case 1471: return "129I";
		case 1472: return "130I";
		case 1473: return "131I";
		case 1474: return "132I";
		case 1475: return "133I";
		case 1476: return "134I";
		case 1477: return "135I";
		case 1478: return "136I";
		case 1479: return "137I";
		case 1480: return "138I";
		case 1481: return "139I";
		case 1482: return "140I";
		case 1483: return "141I";
		case 1484: return "142I";
		case 1485: return "143I";
		case 1486: return "144I";
		case 1487: return "145I";
		case 1488: return "109Xe";
		case 1489: return "110Xe";
		case 1490: return "111Xe";
		case 1491: return "112Xe";
		case 1492: return "113Xe";
		case 1493: return "114Xe";
		case 1494: return "115Xe";
		case 1495: return "116Xe";
		case 1496: return "117Xe";
		case 1497: return "118Xe";
		case 1498: return "119Xe";
		case 1499: return "120Xe";
		case 1500: return "121Xe";
		case 1501: return "122Xe";
		case 1502: return "123Xe";
		case 1503: return "124Xe";
		case 1504: return "125Xe";
		case 1505: return "126Xe";
		case 1506: return "127Xe";
		case 1507: return "128Xe";
		case 1508: return "129Xe";
		case 1509: return "130Xe";
		case 1510: return "131Xe";
		case 1511: return "132Xe";
		case 1512: return "133Xe";
		case 1513: return "134Xe";
		case 1514: return "135Xe";
		case 1515: return "136Xe";
		case 1516: return "137Xe";
		case 1517: return "138Xe";
		case 1518: return "139Xe";
		case 1519: return "140Xe";
		case 1520: return "141Xe";
		case 1521: return "142Xe";
		case 1522: return "143Xe";
		case 1523: return "144Xe";
		case 1524: return "145Xe";
		case 1525: return "146Xe";
		case 1526: return "147Xe";
		case 1527: return "148Xe";
		case 1528: return "112Cs";
		case 1529: return "113Cs";
		case 1530: return "114Cs";
		case 1531: return "115Cs";
		case 1532: return "116Cs";
		case 1533: return "117Cs";
		case 1534: return "118Cs";
		case 1535: return "119Cs";
		case 1536: return "120Cs";
		case 1537: return "121Cs";
		case 1538: return "122Cs";
		case 1539: return "123Cs";
		case 1540: return "124Cs";
		case 1541: return "125Cs";
		case 1542: return "126Cs";
		case 1543: return "127Cs";
		case 1544: return "128Cs";
		case 1545: return "129Cs";
		case 1546: return "130Cs";
		case 1547: return "131Cs";
		case 1548: return "132Cs";
		case 1549: return "133Cs";
		case 1550: return "134Cs";
		case 1551: return "135Cs";
		case 1552: return "136Cs";
		case 1553: return "137Cs";
		case 1554: return "138Cs";
		case 1555: return "139Cs";
		case 1556: return "140Cs";
		case 1557: return "141Cs";
		case 1558: return "142Cs";
		case 1559: return "143Cs";
		case 1560: return "144Cs";
		case 1561: return "145Cs";
		case 1562: return "146Cs";
		case 1563: return "147Cs";
		case 1564: return "148Cs";
		case 1565: return "149Cs";
		case 1566: return "150Cs";
		case 1567: return "151Cs";
		case 1568: return "114Ba";
		case 1569: return "115Ba";
		case 1570: return "116Ba";
		case 1571: return "117Ba";
		case 1572: return "118Ba";
		case 1573: return "119Ba";
		case 1574: return "120Ba";
		case 1575: return "121Ba";
		case 1576: return "122Ba";
		case 1577: return "123Ba";
		case 1578: return "124Ba";
		case 1579: return "125Ba";
		case 1580: return "126Ba";
		case 1581: return "127Ba";
		case 1582: return "128Ba";
		case 1583: return "129Ba";
		case 1584: return "130Ba";
		case 1585: return "131Ba";
		case 1586: return "132Ba";
		case 1587: return "133Ba";
		case 1588: return "134Ba";
		case 1589: return "135Ba";
		case 1590: return "136Ba";
		case 1591: return "137Ba";
		case 1592: return "138Ba";
		case 1593: return "139Ba";
		case 1594: return "140Ba";
		case 1595: return "141Ba";
		case 1596: return "142Ba";
		case 1597: return "143Ba";
		case 1598: return "144Ba";
		case 1599: return "145Ba";
		case 1600: return "146Ba";
		case 1601: return "147Ba";
		case 1602: return "148Ba";
		case 1603: return "149Ba";
		case 1604: return "150Ba";
		case 1605: return "151Ba";
		case 1606: return "152Ba";
		case 1607: return "153Ba";
		case 1608: return "116La";
		case 1609: return "117La";
		case 1610: return "118La";
		case 1611: return "119La";
		case 1612: return "120La";
		case 1613: return "121La";
		case 1614: return "122La";
		case 1615: return "123La";
		case 1616: return "124La";
		case 1617: return "125La";
		case 1618: return "126La";
		case 1619: return "127La";
		case 1620: return "128La";
		case 1621: return "129La";
		case 1622: return "130La";
		case 1623: return "131La";
		case 1624: return "132La";
		case 1625: return "133La";
		case 1626: return "134La";
		case 1627: return "135La";
		case 1628: return "136La";
		case 1629: return "137La";
		case 1630: return "138La";
		case 1631: return "139La";
		case 1632: return "140La";
		case 1633: return "141La";
		case 1634: return "142La";
		case 1635: return "143La";
		case 1636: return "144La";
		case 1637: return "145La";
		case 1638: return "146La";
		case 1639: return "147La";
		case 1640: return "148La";
		case 1641: return "149La";
		case 1642: return "150La";
		case 1643: return "151La";
		case 1644: return "152La";
		case 1645: return "153La";
		case 1646: return "154La";
		case 1647: return "155La";
		case 1648: return "119Ce";
		case 1649: return "120Ce";
		case 1650: return "121Ce";
		case 1651: return "122Ce";
		case 1652: return "123Ce";
		case 1653: return "124Ce";
		case 1654: return "125Ce";
		case 1655: return "126Ce";
		case 1656: return "127Ce";
		case 1657: return "128Ce";
		case 1658: return "129Ce";
		case 1659: return "130Ce";
		case 1660: return "131Ce";
		case 1661: return "132Ce";
		case 1662: return "133Ce";
		case 1663: return "134Ce";
		case 1664: return "135Ce";
		case 1665: return "136Ce";
		case 1666: return "137Ce";
		case 1667: return "138Ce";
		case 1668: return "139Ce";
		case 1669: return "140Ce";
		case 1670: return "141Ce";
		case 1671: return "142Ce";
		case 1672: return "143Ce";
		case 1673: return "144Ce";
		case 1674: return "145Ce";
		case 1675: return "146Ce";
		case 1676: return "147Ce";
		case 1677: return "148Ce";
		case 1678: return "149Ce";
		case 1679: return "150Ce";
		case 1680: return "151Ce";
		case 1681: return "152Ce";
		case 1682: return "153Ce";
		case 1683: return "154Ce";
		case 1684: return "155Ce";
		case 1685: return "156Ce";
		case 1686: return "157Ce";
		case 1687: return "121Pr";
		case 1688: return "122Pr";
		case 1689: return "123Pr";
		case 1690: return "124Pr";
		case 1691: return "125Pr";
		case 1692: return "126Pr";
		case 1693: return "127Pr";
		case 1694: return "128Pr";
		case 1695: return "129Pr";
		case 1696: return "130Pr";
		case 1697: return "131Pr";
		case 1698: return "132Pr";
		case 1699: return "133Pr";
		case 1700: return "134Pr";
		case 1701: return "135Pr";
		case 1702: return "136Pr";
		case 1703: return "137Pr";
		case 1704: return "138Pr";
		case 1705: return "139Pr";
		case 1706: return "140Pr";
		case 1707: return "141Pr";
		case 1708: return "142Pr";
		case 1709: return "143Pr";
		case 1710: return "144Pr";
		case 1711: return "145Pr";
		case 1712: return "146Pr";
		case 1713: return "147Pr";
		case 1714: return "148Pr";
		case 1715: return "149Pr";
		case 1716: return "150Pr";
		case 1717: return "151Pr";
		case 1718: return "152Pr";
		case 1719: return "153Pr";
		case 1720: return "154Pr";
		case 1721: return "155Pr";
		case 1722: return "156Pr";
		case 1723: return "157Pr";
		case 1724: return "158Pr";
		case 1725: return "159Pr";
		case 1726: return "124Nd";
		case 1727: return "125Nd";
		case 1728: return "126Nd";
		case 1729: return "127Nd";
		case 1730: return "128Nd";
		case 1731: return "129Nd";
		case 1732: return "130Nd";
		case 1733: return "131Nd";
		case 1734: return "132Nd";
		case 1735: return "133Nd";
		case 1736: return "134Nd";
		case 1737: return "135Nd";
		case 1738: return "136Nd";
		case 1739: return "137Nd";
		case 1740: return "138Nd";
		case 1741: return "139Nd";
		case 1742: return "140Nd";
		case 1743: return "141Nd";
		case 1744: return "142Nd";
		case 1745: return "143Nd";
		case 1746: return "144Nd";
		case 1747: return "145Nd";
		case 1748: return "146Nd";
		case 1749: return "147Nd";
		case 1750: return "148Nd";
		case 1751: return "149Nd";
		case 1752: return "150Nd";
		case 1753: return "151Nd";
		case 1754: return "152Nd";
		case 1755: return "153Nd";
		case 1756: return "154Nd";
		case 1757: return "155Nd";
		case 1758: return "156Nd";
		case 1759: return "157Nd";
		case 1760: return "158Nd";
		case 1761: return "159Nd";
		case 1762: return "160Nd";
		case 1763: return "161Nd";
		case 1764: return "126Pm";
		case 1765: return "127Pm";
		case 1766: return "128Pm";
		case 1767: return "129Pm";
		case 1768: return "130Pm";
		case 1769: return "131Pm";
		case 1770: return "132Pm";
		case 1771: return "133Pm";
		case 1772: return "134Pm";
		case 1773: return "135Pm";
		case 1774: return "136Pm";
		case 1775: return "137Pm";
		case 1776: return "138Pm";
		case 1777: return "139Pm";
		case 1778: return "140Pm";
		case 1779: return "141Pm";
		case 1780: return "142Pm";
		case 1781: return "143Pm";
		case 1782: return "144Pm";
		case 1783: return "145Pm";
		case 1784: return "146Pm";
		case 1785: return "147Pm";
		case 1786: return "148Pm";
		case 1787: return "149Pm";
		case 1788: return "150Pm";
		case 1789: return "151Pm";
		case 1790: return "152Pm";
		case 1791: return "153Pm";
		case 1792: return "154Pm";
		case 1793: return "155Pm";
		case 1794: return "156Pm";
		case 1795: return "157Pm";
		case 1796: return "158Pm";
		case 1797: return "159Pm";
		case 1798: return "160Pm";
		case 1799: return "161Pm";
		case 1800: return "162Pm";
		case 1801: return "163Pm";
		case 1802: return "128Sm";
		case 1803: return "129Sm";
		case 1804: return "130Sm";
		case 1805: return "131Sm";
		case 1806: return "132Sm";
		case 1807: return "133Sm";
		case 1808: return "134Sm";
		case 1809: return "135Sm";
		case 1810: return "136Sm";
		case 1811: return "137Sm";
		case 1812: return "138Sm";
		case 1813: return "139Sm";
		case 1814: return "140Sm";
		case 1815: return "141Sm";
		case 1816: return "142Sm";
		case 1817: return "143Sm";
		case 1818: return "144Sm";
		case 1819: return "145Sm";
		case 1820: return "146Sm";
		case 1821: return "147Sm";
		case 1822: return "148Sm";
		case 1823: return "149Sm";
		case 1824: return "150Sm";
		case 1825: return "151Sm";
		case 1826: return "152Sm";
		case 1827: return "153Sm";
		case 1828: return "154Sm";
		case 1829: return "155Sm";
		case 1830: return "156Sm";
		case 1831: return "157Sm";
		case 1832: return "158Sm";
		case 1833: return "159Sm";
		case 1834: return "160Sm";
		case 1835: return "161Sm";
		case 1836: return "162Sm";
		case 1837: return "163Sm";
		case 1838: return "164Sm";
		case 1839: return "165Sm";
		case 1840: return "130Eu";
		case 1841: return "131Eu";
		case 1842: return "132Eu";
		case 1843: return "133Eu";
		case 1844: return "134Eu";
		case 1845: return "135Eu";
		case 1846: return "136Eu";
		case 1847: return "137Eu";
		case 1848: return "138Eu";
		case 1849: return "139Eu";
		case 1850: return "140Eu";
		case 1851: return "141Eu";
		case 1852: return "142Eu";
		case 1853: return "143Eu";
		case 1854: return "144Eu";
		case 1855: return "145Eu";
		case 1856: return "146Eu";
		case 1857: return "147Eu";
		case 1858: return "148Eu";
		case 1859: return "149Eu";
		case 1860: return "150Eu";
		case 1861: return "151Eu";
		case 1862: return "152Eu";
		case 1863: return "153Eu";
		case 1864: return "154Eu";
		case 1865: return "155Eu";
		case 1866: return "156Eu";
		case 1867: return "157Eu";
		case 1868: return "158Eu";
		case 1869: return "159Eu";
		case 1870: return "160Eu";
		case 1871: return "161Eu";
		case 1872: return "162Eu";
		case 1873: return "163Eu";
		case 1874: return "164Eu";
		case 1875: return "165Eu";
		case 1876: return "166Eu";
		case 1877: return "167Eu";
		case 1878: return "133Gd";
		case 1879: return "134Gd";
		case 1880: return "135Gd";
		case 1881: return "136Gd";
		case 1882: return "137Gd";
		case 1883: return "138Gd";
		case 1884: return "139Gd";
		case 1885: return "140Gd";
		case 1886: return "141Gd";
		case 1887: return "142Gd";
		case 1888: return "143Gd";
		case 1889: return "144Gd";
		case 1890: return "145Gd";
		case 1891: return "146Gd";
		case 1892: return "147Gd";
		case 1893: return "148Gd";
		case 1894: return "149Gd";
		case 1895: return "150Gd";
		case 1896: return "151Gd";
		case 1897: return "152Gd";
		case 1898: return "153Gd";
		case 1899: return "154Gd";
		case 1900: return "155Gd";
		case 1901: return "156Gd";
		case 1902: return "157Gd";
		case 1903: return "158Gd";
		case 1904: return "159Gd";
		case 1905: return "160Gd";
		case 1906: return "161Gd";
		case 1907: return "162Gd";
		case 1908: return "163Gd";
		case 1909: return "164Gd";
		case 1910: return "165Gd";
		case 1911: return "166Gd";
		case 1912: return "167Gd";
		case 1913: return "168Gd";
		case 1914: return "169Gd";
		case 1915: return "135Tb";
		case 1916: return "136Tb";
		case 1917: return "137Tb";
		case 1918: return "138Tb";
		case 1919: return "139Tb";
		case 1920: return "140Tb";
		case 1921: return "141Tb";
		case 1922: return "142Tb";
		case 1923: return "143Tb";
		case 1924: return "144Tb";
		case 1925: return "145Tb";
		case 1926: return "146Tb";
		case 1927: return "147Tb";
		case 1928: return "148Tb";
		case 1929: return "149Tb";
		case 1930: return "150Tb";
		case 1931: return "151Tb";
		case 1932: return "152Tb";
		case 1933: return "153Tb";
		case 1934: return "154Tb";
		case 1935: return "155Tb";
		case 1936: return "156Tb";
		case 1937: return "157Tb";
		case 1938: return "158Tb";
		case 1939: return "159Tb";
		case 1940: return "160Tb";
		case 1941: return "161Tb";
		case 1942: return "162Tb";
		case 1943: return "163Tb";
		case 1944: return "164Tb";
		case 1945: return "165Tb";
		case 1946: return "166Tb";
		case 1947: return "167Tb";
		case 1948: return "168Tb";
		case 1949: return "169Tb";
		case 1950: return "170Tb";
		case 1951: return "171Tb";
		case 1952: return "138Dy";
		case 1953: return "139Dy";
		case 1954: return "140Dy";
		case 1955: return "141Dy";
		case 1956: return "142Dy";
		case 1957: return "143Dy";
		case 1958: return "144Dy";
		case 1959: return "145Dy";
		case 1960: return "146Dy";
		case 1961: return "147Dy";
		case 1962: return "148Dy";
		case 1963: return "149Dy";
		case 1964: return "150Dy";
		case 1965: return "151Dy";
		case 1966: return "152Dy";
		case 1967: return "153Dy";
		case 1968: return "154Dy";
		case 1969: return "155Dy";
		case 1970: return "156Dy";
		case 1971: return "157Dy";
		case 1972: return "158Dy";
		case 1973: return "159Dy";
		case 1974: return "160Dy";
		case 1975: return "161Dy";
		case 1976: return "162Dy";
		case 1977: return "163Dy";
		case 1978: return "164Dy";
		case 1979: return "165Dy";
		case 1980: return "166Dy";
		case 1981: return "167Dy";
		case 1982: return "168Dy";
		case 1983: return "169Dy";
		case 1984: return "170Dy";
		case 1985: return "171Dy";
		case 1986: return "172Dy";
		case 1987: return "173Dy";
		case 1988: return "140Ho";
		case 1989: return "141Ho";
		case 1990: return "142Ho";
		case 1991: return "143Ho";
		case 1992: return "144Ho";
		case 1993: return "145Ho";
		case 1994: return "146Ho";
		case 1995: return "147Ho";
		case 1996: return "148Ho";
		case 1997: return "149Ho";
		case 1998: return "150Ho";
		case 1999: return "151Ho";
		case 2000: return "152Ho";
		case 2001: return "153Ho";
		case 2002: return "154Ho";
		case 2003: return "155Ho";
		case 2004: return "156Ho";
		case 2005: return "157Ho";
		case 2006: return "158Ho";
		case 2007: return "159Ho";
		case 2008: return "160Ho";
		case 2009: return "161Ho";
		case 2010: return "162Ho";
		case 2011: return "163Ho";
		case 2012: return "164Ho";
		case 2013: return "165Ho";
		case 2014: return "166Ho";
		case 2015: return "167Ho";
		case 2016: return "168Ho";
		case 2017: return "169Ho";
		case 2018: return "170Ho";
		case 2019: return "171Ho";
		case 2020: return "172Ho";
		case 2021: return "173Ho";
		case 2022: return "174Ho";
		case 2023: return "175Ho";
		case 2024: return "142Er";
		case 2025: return "143Er";
		case 2026: return "144Er";
		case 2027: return "145Er";
		case 2028: return "146Er";
		case 2029: return "147Er";
		case 2030: return "148Er";
		case 2031: return "149Er";
		case 2032: return "150Er";
		case 2033: return "151Er";
		case 2034: return "152Er";
		case 2035: return "153Er";
		case 2036: return "154Er";
		case 2037: return "155Er";
		case 2038: return "156Er";
		case 2039: return "157Er";
		case 2040: return "158Er";
		case 2041: return "159Er";
		case 2042: return "160Er";
		case 2043: return "161Er";
		case 2044: return "162Er";
		case 2045: return "163Er";
		case 2046: return "164Er";
		case 2047: return "165Er";
		case 2048: return "166Er";
		case 2049: return "167Er";
		case 2050: return "168Er";
		case 2051: return "169Er";
		case 2052: return "170Er";
		case 2053: return "171Er";
		case 2054: return "172Er";
		case 2055: return "173Er";
		case 2056: return "174Er";
		case 2057: return "175Er";
		case 2058: return "176Er";
		case 2059: return "177Er";
		case 2060: return "144Tm";
		case 2061: return "145Tm";
		case 2062: return "146Tm";
		case 2063: return "147Tm";
		case 2064: return "148Tm";
		case 2065: return "149Tm";
		case 2066: return "150Tm";
		case 2067: return "151Tm";
		case 2068: return "152Tm";
		case 2069: return "153Tm";
		case 2070: return "154Tm";
		case 2071: return "155Tm";
		case 2072: return "156Tm";
		case 2073: return "157Tm";
		case 2074: return "158Tm";
		case 2075: return "159Tm";
		case 2076: return "160Tm";
		case 2077: return "161Tm";
		case 2078: return "162Tm";
		case 2079: return "163Tm";
		case 2080: return "164Tm";
		case 2081: return "165Tm";
		case 2082: return "166Tm";
		case 2083: return "167Tm";
		case 2084: return "168Tm";
		case 2085: return "169Tm";
		case 2086: return "170Tm";
		case 2087: return "171Tm";
		case 2088: return "172Tm";
		case 2089: return "173Tm";
		case 2090: return "174Tm";
		case 2091: return "175Tm";
		case 2092: return "176Tm";
		case 2093: return "177Tm";
		case 2094: return "178Tm";
		case 2095: return "179Tm";
		case 2096: return "148Yb";
		case 2097: return "149Yb";
		case 2098: return "150Yb";
		case 2099: return "151Yb";
		case 2100: return "152Yb";
		case 2101: return "153Yb";
		case 2102: return "154Yb";
		case 2103: return "155Yb";
		case 2104: return "156Yb";
		case 2105: return "157Yb";
		case 2106: return "158Yb";
		case 2107: return "159Yb";
		case 2108: return "160Yb";
		case 2109: return "161Yb";
		case 2110: return "162Yb";
		case 2111: return "163Yb";
		case 2112: return "164Yb";
		case 2113: return "165Yb";
		case 2114: return "166Yb";
		case 2115: return "167Yb";
		case 2116: return "168Yb";
		case 2117: return "169Yb";
		case 2118: return "170Yb";
		case 2119: return "171Yb";
		case 2120: return "172Yb";
		case 2121: return "173Yb";
		case 2122: return "174Yb";
		case 2123: return "175Yb";
		case 2124: return "176Yb";
		case 2125: return "177Yb";
		case 2126: return "178Yb";
		case 2127: return "179Yb";
		case 2128: return "180Yb";
		case 2129: return "181Yb";
		case 2130: return "150Lu";
		case 2131: return "151Lu";
		case 2132: return "152Lu";
		case 2133: return "153Lu";
		case 2134: return "154Lu";
		case 2135: return "155Lu";
		case 2136: return "156Lu";
		case 2137: return "157Lu";
		case 2138: return "158Lu";
		case 2139: return "159Lu";
		case 2140: return "160Lu";
		case 2141: return "161Lu";
		case 2142: return "162Lu";
		case 2143: return "163Lu";
		case 2144: return "164Lu";
		case 2145: return "165Lu";
		case 2146: return "166Lu";
		case 2147: return "167Lu";
		case 2148: return "168Lu";
		case 2149: return "169Lu";
		case 2150: return "170Lu";
		case 2151: return "171Lu";
		case 2152: return "172Lu";
		case 2153: return "173Lu";
		case 2154: return "174Lu";
		case 2155: return "175Lu";
		case 2156: return "176Lu";
		case 2157: return "177Lu";
		case 2158: return "178Lu";
		case 2159: return "179Lu";
		case 2160: return "180Lu";
		case 2161: return "181Lu";
		case 2162: return "182Lu";
		case 2163: return "183Lu";
		case 2164: return "184Lu";
		case 2165: return "185Lu";
		case 2166: return "153Hf";
		case 2167: return "154Hf";
		case 2168: return "155Hf";
		case 2169: return "156Hf";
		case 2170: return "157Hf";
		case 2171: return "158Hf";
		case 2172: return "159Hf";
		case 2173: return "160Hf";
		case 2174: return "161Hf";
		case 2175: return "162Hf";
		case 2176: return "163Hf";
		case 2177: return "164Hf";
		case 2178: return "165Hf";
		case 2179: return "166Hf";
		case 2180: return "167Hf";
		case 2181: return "168Hf";
		case 2182: return "169Hf";
		case 2183: return "170Hf";
		case 2184: return "171Hf";
		case 2185: return "172Hf";
		case 2186: return "173Hf";
		case 2187: return "174Hf";
		case 2188: return "175Hf";
		case 2189: return "176Hf";
		case 2190: return "177Hf";
		case 2191: return "178Hf";
		case 2192: return "179Hf";
		case 2193: return "180Hf";
		case 2194: return "181Hf";
		case 2195: return "182Hf";
		case 2196: return "183Hf";
		case 2197: return "184Hf";
		case 2198: return "185Hf";
		case 2199: return "186Hf";
		case 2200: return "187Hf";
		case 2201: return "188Hf";
		case 2202: return "189Hf";
		case 2203: return "155Ta";
		case 2204: return "156Ta";
		case 2205: return "157Ta";
		case 2206: return "158Ta";
		case 2207: return "159Ta";
		case 2208: return "160Ta";
		case 2209: return "161Ta";
		case 2210: return "162Ta";
		case 2211: return "163Ta";
		case 2212: return "164Ta";
		case 2213: return "165Ta";
		case 2214: return "166Ta";
		case 2215: return "167Ta";
		case 2216: return "168Ta";
		case 2217: return "169Ta";
		case 2218: return "170Ta";
		case 2219: return "171Ta";
		case 2220: return "172Ta";
		case 2221: return "173Ta";
		case 2222: return "174Ta";
		case 2223: return "175Ta";
		case 2224: return "176Ta";
		case 2225: return "177Ta";
		case 2226: return "178Ta";
		case 2227: return "179Ta";
		case 2228: return "180Ta";
		case 2229: return "181Ta";
		case 2230: return "182Ta";
		case 2231: return "183Ta";
		case 2232: return "184Ta";
		case 2233: return "185Ta";
		case 2234: return "186Ta";
		case 2235: return "187Ta";
		case 2236: return "188Ta";
		case 2237: return "189Ta";
		case 2238: return "190Ta";
		case 2239: return "191Ta";
		case 2240: return "192Ta";
		case 2241: return "157W";
		case 2242: return "158W";
		case 2243: return "159W";
		case 2244: return "160W";
		case 2245: return "161W";
		case 2246: return "162W";
		case 2247: return "163W";
		case 2248: return "164W";
		case 2249: return "165W";
		case 2250: return "166W";
		case 2251: return "167W";
		case 2252: return "168W";
		case 2253: return "169W";
		case 2254: return "170W";
		case 2255: return "171W";
		case 2256: return "172W";
		case 2257: return "173W";
		case 2258: return "174W";
		case 2259: return "175W";
		case 2260: return "176W";
		case 2261: return "177W";
		case 2262: return "178W";
		case 2263: return "179W";
		case 2264: return "180W";
		case 2265: return "181W";
		case 2266: return "182W";
		case 2267: return "183W";
		case 2268: return "184W";
		case 2269: return "185W";
		case 2270: return "186W";
		case 2271: return "187W";
		case 2272: return "188W";
		case 2273: return "189W";
		case 2274: return "190W";
		case 2275: return "191W";
		case 2276: return "192W";
		case 2277: return "193W";
		case 2278: return "194W";
		case 2279: return "159Re";
		case 2280: return "160Re";
		case 2281: return "161Re";
		case 2282: return "162Re";
		case 2283: return "163Re";
		case 2284: return "164Re";
		case 2285: return "165Re";
		case 2286: return "166Re";
		case 2287: return "167Re";
		case 2288: return "168Re";
		case 2289: return "169Re";
		case 2290: return "170Re";
		case 2291: return "171Re";
		case 2292: return "172Re";
		case 2293: return "173Re";
		case 2294: return "174Re";
		case 2295: return "175Re";
		case 2296: return "176Re";
		case 2297: return "177Re";
		case 2298: return "178Re";
		case 2299: return "179Re";
		case 2300: return "180Re";
		case 2301: return "181Re";
		case 2302: return "182Re";
		case 2303: return "183Re";
		case 2304: return "184Re";
		case 2305: return "185Re";
		case 2306: return "186Re";
		case 2307: return "187Re";
		case 2308: return "188Re";
		case 2309: return "189Re";
		case 2310: return "190Re";
		case 2311: return "191Re";
		case 2312: return "192Re";
		case 2313: return "193Re";
		case 2314: return "194Re";
		case 2315: return "195Re";
		case 2316: return "196Re";
		case 2317: return "197Re";
		case 2318: return "198Re";
		case 2319: return "161Os";
		case 2320: return "162Os";
		case 2321: return "163Os";
		case 2322: return "164Os";
		case 2323: return "165Os";
		case 2324: return "166Os";
		case 2325: return "167Os";
		case 2326: return "168Os";
		case 2327: return "169Os";
		case 2328: return "170Os";
		case 2329: return "171Os";
		case 2330: return "172Os";
		case 2331: return "173Os";
		case 2332: return "174Os";
		case 2333: return "175Os";
		case 2334: return "176Os";
		case 2335: return "177Os";
		case 2336: return "178Os";
		case 2337: return "179Os";
		case 2338: return "180Os";
		case 2339: return "181Os";
		case 2340: return "182Os";
		case 2341: return "183Os";
		case 2342: return "184Os";
		case 2343: return "185Os";
		case 2344: return "186Os";
		case 2345: return "187Os";
		case 2346: return "188Os";
		case 2347: return "189Os";
		case 2348: return "190Os";
		case 2349: return "191Os";
		case 2350: return "192Os";
		case 2351: return "193Os";
		case 2352: return "194Os";
		case 2353: return "195Os";
		case 2354: return "196Os";
		case 2355: return "197Os";
		case 2356: return "198Os";
		case 2357: return "199Os";
		case 2358: return "200Os";
		case 2359: return "201Os";
		case 2360: return "202Os";
		case 2361: return "164Ir";
		case 2362: return "165Ir";
		case 2363: return "166Ir";
		case 2364: return "167Ir";
		case 2365: return "168Ir";
		case 2366: return "169Ir";
		case 2367: return "170Ir";
		case 2368: return "171Ir";
		case 2369: return "172Ir";
		case 2370: return "173Ir";
		case 2371: return "174Ir";
		case 2372: return "175Ir";
		case 2373: return "176Ir";
		case 2374: return "177Ir";
		case 2375: return "178Ir";
		case 2376: return "179Ir";
		case 2377: return "180Ir";
		case 2378: return "181Ir";
		case 2379: return "182Ir";
		case 2380: return "183Ir";
		case 2381: return "184Ir";
		case 2382: return "185Ir";
		case 2383: return "186Ir";
		case 2384: return "187Ir";
		case 2385: return "188Ir";
		case 2386: return "189Ir";
		case 2387: return "190Ir";
		case 2388: return "191Ir";
		case 2389: return "192Ir";
		case 2390: return "193Ir";
		case 2391: return "194Ir";
		case 2392: return "195Ir";
		case 2393: return "196Ir";
		case 2394: return "197Ir";
		case 2395: return "198Ir";
		case 2396: return "199Ir";
		case 2397: return "200Ir";
		case 2398: return "201Ir";
		case 2399: return "202Ir";
		case 2400: return "203Ir";
		case 2401: return "204Ir";
		case 2402: return "166Pt";
		case 2403: return "167Pt";
		case 2404: return "168Pt";
		case 2405: return "169Pt";
		case 2406: return "170Pt";
		case 2407: return "171Pt";
		case 2408: return "172Pt";
		case 2409: return "173Pt";
		case 2410: return "174Pt";
		case 2411: return "175Pt";
		case 2412: return "176Pt";
		case 2413: return "177Pt";
		case 2414: return "178Pt";
		case 2415: return "179Pt";
		case 2416: return "180Pt";
		case 2417: return "181Pt";
		case 2418: return "182Pt";
		case 2419: return "183Pt";
		case 2420: return "184Pt";
		case 2421: return "185Pt";
		case 2422: return "186Pt";
		case 2423: return "187Pt";
		case 2424: return "188Pt";
		case 2425: return "189Pt";
		case 2426: return "190Pt";
		case 2427: return "191Pt";
		case 2428: return "192Pt";
		case 2429: return "193Pt";
		case 2430: return "194Pt";
		case 2431: return "195Pt";
		case 2432: return "196Pt";
		case 2433: return "197Pt";
		case 2434: return "198Pt";
		case 2435: return "199Pt";
		case 2436: return "200Pt";
		case 2437: return "201Pt";
		case 2438: return "202Pt";
		case 2439: return "203Pt";
		case 2440: return "204Pt";
		case 2441: return "205Pt";
		case 2442: return "206Pt";
		case 2443: return "169Au";
		case 2444: return "170Au";
		case 2445: return "171Au";
		case 2446: return "172Au";
		case 2447: return "173Au";
		case 2448: return "174Au";
		case 2449: return "175Au";
		case 2450: return "176Au";
		case 2451: return "177Au";
		case 2452: return "178Au";
		case 2453: return "179Au";
		case 2454: return "180Au";
		case 2455: return "181Au";
		case 2456: return "182Au";
		case 2457: return "183Au";
		case 2458: return "184Au";
		case 2459: return "185Au";
		case 2460: return "186Au";
		case 2461: return "187Au";
		case 2462: return "188Au";
		case 2463: return "189Au";
		case 2464: return "190Au";
		case 2465: return "191Au";
		case 2466: return "192Au";
		case 2467: return "193Au";
		case 2468: return "194Au";
		case 2469: return "195Au";
		case 2470: return "196Au";
		case 2471: return "197Au";
		case 2472: return "198Au";
		case 2473: return "199Au";
		case 2474: return "200Au";
		case 2475: return "201Au";
		case 2476: return "202Au";
		case 2477: return "203Au";
		case 2478: return "204Au";
		case 2479: return "205Au";
		case 2480: return "206Au";
		case 2481: return "207Au";
		case 2482: return "208Au";
		case 2483: return "209Au";
		case 2484: return "210Au";
		case 2485: return "171Hg";
		case 2486: return "172Hg";
		case 2487: return "173Hg";
		case 2488: return "174Hg";
		case 2489: return "175Hg";
		case 2490: return "176Hg";
		case 2491: return "177Hg";
		case 2492: return "178Hg";
		case 2493: return "179Hg";
		case 2494: return "180Hg";
		case 2495: return "181Hg";
		case 2496: return "182Hg";
		case 2497: return "183Hg";
		case 2498: return "184Hg";
		case 2499: return "185Hg";
		case 2500: return "186Hg";
		case 2501: return "187Hg";
		case 2502: return "188Hg";
		case 2503: return "189Hg";
		case 2504: return "190Hg";
		case 2505: return "191Hg";
		case 2506: return "192Hg";
		case 2507: return "193Hg";
		case 2508: return "194Hg";
		case 2509: return "195Hg";
		case 2510: return "196Hg";
		case 2511: return "197Hg";
		case 2512: return "198Hg";
		case 2513: return "199Hg";
		case 2514: return "200Hg";
		case 2515: return "201Hg";
		case 2516: return "202Hg";
		case 2517: return "203Hg";
		case 2518: return "204Hg";
		case 2519: return "205Hg";
		case 2520: return "206Hg";
		case 2521: return "207Hg";
		case 2522: return "208Hg";
		case 2523: return "209Hg";
		case 2524: return "210Hg";
		case 2525: return "211Hg";
		case 2526: return "212Hg";
		case 2527: return "213Hg";
		case 2528: return "214Hg";
		case 2529: return "215Hg";
		case 2530: return "216Hg";
		case 2531: return "176Tl";
		case 2532: return "177Tl";
		case 2533: return "178Tl";
		case 2534: return "179Tl";
		case 2535: return "180Tl";
		case 2536: return "181Tl";
		case 2537: return "182Tl";
		case 2538: return "183Tl";
		case 2539: return "184Tl";
		case 2540: return "185Tl";
		case 2541: return "186Tl";
		case 2542: return "187Tl";
		case 2543: return "188Tl";
		case 2544: return "189Tl";
		case 2545: return "190Tl";
		case 2546: return "191Tl";
		case 2547: return "192Tl";
		case 2548: return "193Tl";
		case 2549: return "194Tl";
		case 2550: return "195Tl";
		case 2551: return "196Tl";
		case 2552: return "197Tl";
		case 2553: return "198Tl";
		case 2554: return "199Tl";
		case 2555: return "200Tl";
		case 2556: return "201Tl";
		case 2557: return "202Tl";
		case 2558: return "203Tl";
		case 2559: return "204Tl";
		case 2560: return "205Tl";
		case 2561: return "206Tl";
		case 2562: return "207Tl";
		case 2563: return "208Tl";
		case 2564: return "209Tl";
		case 2565: return "210Tl";
		case 2566: return "211Tl";
		case 2567: return "212Tl";
		case 2568: return "213Tl";
		case 2569: return "214Tl";
		case 2570: return "215Tl";
		case 2571: return "216Tl";
		case 2572: return "217Tl";
		case 2573: return "218Tl";
		case 2574: return "178Pb";
		case 2575: return "179Pb";
		case 2576: return "180Pb";
		case 2577: return "181Pb";
		case 2578: return "182Pb";
		case 2579: return "183Pb";
		case 2580: return "184Pb";
		case 2581: return "185Pb";
		case 2582: return "186Pb";
		case 2583: return "187Pb";
		case 2584: return "188Pb";
		case 2585: return "189Pb";
		case 2586: return "190Pb";
		case 2587: return "191Pb";
		case 2588: return "192Pb";
		case 2589: return "193Pb";
		case 2590: return "194Pb";
		case 2591: return "195Pb";
		case 2592: return "196Pb";
		case 2593: return "197Pb";
		case 2594: return "198Pb";
		case 2595: return "199Pb";
		case 2596: return "200Pb";
		case 2597: return "201Pb";
		case 2598: return "202Pb";
		case 2599: return "203Pb";
		case 2600: return "204Pb";
		case 2601: return "205Pb";
		case 2602: return "206Pb";
		case 2603: return "207Pb";
		case 2604: return "208Pb";
		case 2605: return "209Pb";
		case 2606: return "210Pb";
		case 2607: return "211Pb";
		case 2608: return "212Pb";
		case 2609: return "213Pb";
		case 2610: return "214Pb";
		case 2611: return "215Pb";
		case 2612: return "216Pb";
		case 2613: return "217Pb";
		case 2614: return "218Pb";
		case 2615: return "219Pb";
		case 2616: return "220Pb";
		case 2617: return "184Bi";
		case 2618: return "185Bi";
		case 2619: return "186Bi";
		case 2620: return "187Bi";
		case 2621: return "188Bi";
		case 2622: return "189Bi";
		case 2623: return "190Bi";
		case 2624: return "191Bi";
		case 2625: return "192Bi";
		case 2626: return "193Bi";
		case 2627: return "194Bi";
		case 2628: return "195Bi";
		case 2629: return "196Bi";
		case 2630: return "197Bi";
		case 2631: return "198Bi";
		case 2632: return "199Bi";
		case 2633: return "200Bi";
		case 2634: return "201Bi";
		case 2635: return "202Bi";
		case 2636: return "203Bi";
		case 2637: return "204Bi";
		case 2638: return "205Bi";
		case 2639: return "206Bi";
		case 2640: return "207Bi";
		case 2641: return "208Bi";
		case 2642: return "209Bi";
		case 2643: return "210Bi";
		case 2644: return "211Bi";
		case 2645: return "212Bi";
		case 2646: return "213Bi";
		case 2647: return "214Bi";
		case 2648: return "215Bi";
		case 2649: return "216Bi";
		case 2650: return "217Bi";
		case 2651: return "218Bi";
		case 2652: return "219Bi";
		case 2653: return "220Bi";
		case 2654: return "221Bi";
		case 2655: return "222Bi";
		case 2656: return "223Bi";
		case 2657: return "224Bi";
		case 2658: return "186Po";
		case 2659: return "187Po";
		case 2660: return "188Po";
		case 2661: return "189Po";
		case 2662: return "190Po";
		case 2663: return "191Po";
		case 2664: return "192Po";
		case 2665: return "193Po";
		case 2666: return "194Po";
		case 2667: return "195Po";
		case 2668: return "196Po";
		case 2669: return "197Po";
		case 2670: return "198Po";
		case 2671: return "199Po";
		case 2672: return "200Po";
		case 2673: return "201Po";
		case 2674: return "202Po";
		case 2675: return "203Po";
		case 2676: return "204Po";
		case 2677: return "205Po";
		case 2678: return "206Po";
		case 2679: return "207Po";
		case 2680: return "208Po";
		case 2681: return "209Po";
		case 2682: return "210Po";
		case 2683: return "211Po";
		case 2684: return "212Po";
		case 2685: return "213Po";
		case 2686: return "214Po";
		case 2687: return "215Po";
		case 2688: return "216Po";
		case 2689: return "217Po";
		case 2690: return "218Po";
		case 2691: return "219Po";
		case 2692: return "220Po";
		case 2693: return "221Po";
		case 2694: return "222Po";
		case 2695: return "223Po";
		case 2696: return "224Po";
		case 2697: return "225Po";
		case 2698: return "226Po";
		case 2699: return "227Po";
		case 2700: return "191At";
		case 2701: return "192At";
		case 2702: return "193At";
		case 2703: return "194At";
		case 2704: return "195At";
		case 2705: return "196At";
		case 2706: return "197At";
		case 2707: return "198At";
		case 2708: return "199At";
		case 2709: return "200At";
		case 2710: return "201At";
		case 2711: return "202At";
		case 2712: return "203At";
		case 2713: return "204At";
		case 2714: return "205At";
		case 2715: return "206At";
		case 2716: return "207At";
		case 2717: return "208At";
		case 2718: return "209At";
		case 2719: return "210At";
		case 2720: return "211At";
		case 2721: return "212At";
		case 2722: return "213At";
		case 2723: return "214At";
		case 2724: return "215At";
		case 2725: return "216At";
		case 2726: return "217At";
		case 2727: return "218At";
		case 2728: return "219At";
		case 2729: return "220At";
		case 2730: return "221At";
		case 2731: return "222At";
		case 2732: return "223At";
		case 2733: return "224At";
		case 2734: return "225At";
		case 2735: return "226At";
		case 2736: return "227At";
		case 2737: return "228At";
		case 2738: return "229At";
		case 2739: return "193Rn";
		case 2740: return "194Rn";
		case 2741: return "195Rn";
		case 2742: return "196Rn";
		case 2743: return "197Rn";
		case 2744: return "198Rn";
		case 2745: return "199Rn";
		case 2746: return "200Rn";
		case 2747: return "201Rn";
		case 2748: return "202Rn";
		case 2749: return "203Rn";
		case 2750: return "204Rn";
		case 2751: return "205Rn";
		case 2752: return "206Rn";
		case 2753: return "207Rn";
		case 2754: return "208Rn";
		case 2755: return "209Rn";
		case 2756: return "210Rn";
		case 2757: return "211Rn";
		case 2758: return "212Rn";
		case 2759: return "213Rn";
		case 2760: return "214Rn";
		case 2761: return "215Rn";
		case 2762: return "216Rn";
		case 2763: return "217Rn";
		case 2764: return "218Rn";
		case 2765: return "219Rn";
		case 2766: return "220Rn";
		case 2767: return "221Rn";
		case 2768: return "222Rn";
		case 2769: return "223Rn";
		case 2770: return "224Rn";
		case 2771: return "225Rn";
		case 2772: return "226Rn";
		case 2773: return "227Rn";
		case 2774: return "228Rn";
		case 2775: return "229Rn";
		case 2776: return "230Rn";
		case 2777: return "231Rn";
		case 2778: return "199Fr";
		case 2779: return "200Fr";
		case 2780: return "201Fr";
		case 2781: return "202Fr";
		case 2782: return "203Fr";
		case 2783: return "204Fr";
		case 2784: return "205Fr";
		case 2785: return "206Fr";
		case 2786: return "207Fr";
		case 2787: return "208Fr";
		case 2788: return "209Fr";
		case 2789: return "210Fr";
		case 2790: return "211Fr";
		case 2791: return "212Fr";
		case 2792: return "213Fr";
		case 2793: return "214Fr";
		case 2794: return "215Fr";
		case 2795: return "216Fr";
		case 2796: return "217Fr";
		case 2797: return "218Fr";
		case 2798: return "219Fr";
		case 2799: return "220Fr";
		case 2800: return "221Fr";
		case 2801: return "222Fr";
		case 2802: return "223Fr";
		case 2803: return "224Fr";
		case 2804: return "225Fr";
		case 2805: return "226Fr";
		case 2806: return "227Fr";
		case 2807: return "228Fr";
		case 2808: return "229Fr";
		case 2809: return "230Fr";
		case 2810: return "231Fr";
		case 2811: return "232Fr";
		case 2812: return "233Fr";
		case 2813: return "201Ra";
		case 2814: return "202Ra";
		case 2815: return "203Ra";
		case 2816: return "204Ra";
		case 2817: return "205Ra";
		case 2818: return "206Ra";
		case 2819: return "207Ra";
		case 2820: return "208Ra";
		case 2821: return "209Ra";
		case 2822: return "210Ra";
		case 2823: return "211Ra";
		case 2824: return "212Ra";
		case 2825: return "213Ra";
		case 2826: return "214Ra";
		case 2827: return "215Ra";
		case 2828: return "216Ra";
		case 2829: return "217Ra";
		case 2830: return "218Ra";
		case 2831: return "219Ra";
		case 2832: return "220Ra";
		case 2833: return "221Ra";
		case 2834: return "222Ra";
		case 2835: return "223Ra";
		case 2836: return "224Ra";
		case 2837: return "225Ra";
		case 2838: return "226Ra";
		case 2839: return "227Ra";
		case 2840: return "228Ra";
		case 2841: return "229Ra";
		case 2842: return "230Ra";
		case 2843: return "231Ra";
		case 2844: return "232Ra";
		case 2845: return "233Ra";
		case 2846: return "234Ra";
		case 2847: return "235Ra";
		case 2848: return "206Ac";
		case 2849: return "207Ac";
		case 2850: return "208Ac";
		case 2851: return "209Ac";
		case 2852: return "210Ac";
		case 2853: return "211Ac";
		case 2854: return "212Ac";
		case 2855: return "213Ac";
		case 2856: return "214Ac";
		case 2857: return "215Ac";
		case 2858: return "216Ac";
		case 2859: return "217Ac";
		case 2860: return "218Ac";
		case 2861: return "219Ac";
		case 2862: return "220Ac";
		case 2863: return "221Ac";
		case 2864: return "222Ac";
		case 2865: return "223Ac";
		case 2866: return "224Ac";
		case 2867: return "225Ac";
		case 2868: return "226Ac";
		case 2869: return "227Ac";
		case 2870: return "228Ac";
		case 2871: return "229Ac";
		case 2872: return "230Ac";
		case 2873: return "231Ac";
		case 2874: return "232Ac";
		case 2875: return "233Ac";
		case 2876: return "234Ac";
		case 2877: return "235Ac";
		case 2878: return "236Ac";
		case 2879: return "237Ac";
		case 2880: return "208Th";
		case 2881: return "209Th";
		case 2882: return "210Th";
		case 2883: return "211Th";
		case 2884: return "212Th";
		case 2885: return "213Th";
		case 2886: return "214Th";
		case 2887: return "215Th";
		case 2888: return "216Th";
		case 2889: return "217Th";
		case 2890: return "218Th";
		case 2891: return "219Th";
		case 2892: return "220Th";
		case 2893: return "221Th";
		case 2894: return "222Th";
		case 2895: return "223Th";
		case 2896: return "224Th";
		case 2897: return "225Th";
		case 2898: return "226Th";
		case 2899: return "227Th";
		case 2900: return "228Th";
		case 2901: return "229Th";
		case 2902: return "230Th";
		case 2903: return "231Th";
		case 2904: return "232Th";
		case 2905: return "233Th";
		case 2906: return "234Th";
		case 2907: return "235Th";
		case 2908: return "236Th";
		case 2909: return "237Th";
		case 2910: return "238Th";
		case 2911: return "239Th";
		case 2912: return "212Pa";
		case 2913: return "213Pa";
		case 2914: return "214Pa";
		case 2915: return "215Pa";
		case 2916: return "216Pa";
		case 2917: return "217Pa";
		case 2918: return "218Pa";
		case 2919: return "219Pa";
		case 2920: return "220Pa";
		case 2921: return "221Pa";
		case 2922: return "222Pa";
		case 2923: return "223Pa";
		case 2924: return "224Pa";
		case 2925: return "225Pa";
		case 2926: return "226Pa";
		case 2927: return "227Pa";
		case 2928: return "228Pa";
		case 2929: return "229Pa";
		case 2930: return "230Pa";
		case 2931: return "231Pa";
		case 2932: return "232Pa";
		case 2933: return "233Pa";
		case 2934: return "234Pa";
		case 2935: return "235Pa";
		case 2936: return "236Pa";
		case 2937: return "237Pa";
		case 2938: return "238Pa";
		case 2939: return "239Pa";
		case 2940: return "240Pa";
		case 2941: return "241Pa";
		case 2942: return "217U";
		case 2943: return "218U";
		case 2944: return "219U";
		case 2945: return "220U";
		case 2946: return "221U";
		case 2947: return "222U";
		case 2948: return "223U";
		case 2949: return "224U";
		case 2950: return "225U";
		case 2951: return "226U";
		case 2952: return "227U";
		case 2953: return "228U";
		case 2954: return "229U";
		case 2955: return "230U";
		case 2956: return "231U";
		case 2957: return "232U";
		case 2958: return "233U";
		case 2959: return "234U";
		case 2960: return "235U";
		case 2961: return "236U";
		case 2962: return "237U";
		case 2963: return "238U";
		case 2964: return "239U";
		case 2965: return "240U";
		case 2966: return "241U";
		case 2967: return "242U";
		case 2968: return "243U";
		case 2969: return "219Np";
		case 2970: return "220Np";
		case 2971: return "221Np";
		case 2972: return "222Np";
		case 2973: return "223Np";
		case 2974: return "224Np";
		case 2975: return "225Np";
		case 2976: return "226Np";
		case 2977: return "227Np";
		case 2978: return "228Np";
		case 2979: return "229Np";
		case 2980: return "230Np";
		case 2981: return "231Np";
		case 2982: return "232Np";
		case 2983: return "233Np";
		case 2984: return "234Np";
		case 2985: return "235Np";
		case 2986: return "236Np";
		case 2987: return "237Np";
		case 2988: return "238Np";
		case 2989: return "239Np";
		case 2990: return "240Np";
		case 2991: return "241Np";
		case 2992: return "242Np";
		case 2993: return "243Np";
		case 2994: return "244Np";
		case 2995: return "245Np";
		case 2996: return "228Pu";
		case 2997: return "229Pu";
		case 2998: return "230Pu";
		case 2999: return "231Pu";
		case 3000: return "232Pu";
		case 3001: return "233Pu";
		case 3002: return "234Pu";
		case 3003: return "235Pu";
		case 3004: return "236Pu";
		case 3005: return "237Pu";
		case 3006: return "238Pu";
		case 3007: return "239Pu";
		case 3008: return "240Pu";
		case 3009: return "241Pu";
		case 3010: return "242Pu";
		case 3011: return "243Pu";
		case 3012: return "244Pu";
		case 3013: return "245Pu";
		case 3014: return "246Pu";
		case 3015: return "247Pu";
		case 3016: return "230Am";
		case 3017: return "231Am";
		case 3018: return "232Am";
		case 3019: return "233Am";
		case 3020: return "234Am";
		case 3021: return "235Am";
		case 3022: return "236Am";
		case 3023: return "237Am";
		case 3024: return "238Am";
		case 3025: return "239Am";
		case 3026: return "240Am";
		case 3027: return "241Am";
		case 3028: return "242Am";
		case 3029: return "243Am";
		case 3030: return "244Am";
		case 3031: return "245Am";
		case 3032: return "246Am";
		case 3033: return "247Am";
		case 3034: return "248Am";
		case 3035: return "249Am";
		case 3036: return "232Cm";
		case 3037: return "233Cm";
		case 3038: return "234Cm";
		case 3039: return "235Cm";
		case 3040: return "236Cm";
		case 3041: return "237Cm";
		case 3042: return "238Cm";
		case 3043: return "239Cm";
		case 3044: return "240Cm";
		case 3045: return "241Cm";
		case 3046: return "242Cm";
		case 3047: return "243Cm";
		case 3048: return "244Cm";
		case 3049: return "245Cm";
		case 3050: return "246Cm";
		case 3051: return "247Cm";
		case 3052: return "248Cm";
		case 3053: return "249Cm";
		case 3054: return "250Cm";
		case 3055: return "251Cm";
		case 3056: return "252Cm";
		case 3057: return "234Bk";
		case 3058: return "235Bk";
		case 3059: return "236Bk";
		case 3060: return "237Bk";
		case 3061: return "238Bk";
		case 3062: return "239Bk";
		case 3063: return "240Bk";
		case 3064: return "241Bk";
		case 3065: return "242Bk";
		case 3066: return "243Bk";
		case 3067: return "244Bk";
		case 3068: return "245Bk";
		case 3069: return "246Bk";
		case 3070: return "247Bk";
		case 3071: return "248Bk";
		case 3072: return "249Bk";
		case 3073: return "250Bk";
		case 3074: return "251Bk";
		case 3075: return "252Bk";
		case 3076: return "253Bk";
		case 3077: return "254Bk";
		case 3078: return "237Cf";
		case 3079: return "238Cf";
		case 3080: return "239Cf";
		case 3081: return "240Cf";
		case 3082: return "241Cf";
		case 3083: return "242Cf";
		case 3084: return "243Cf";
		case 3085: return "244Cf";
		case 3086: return "245Cf";
		case 3087: return "246Cf";
		case 3088: return "247Cf";
		case 3089: return "248Cf";
		case 3090: return "249Cf";
		case 3091: return "250Cf";
		case 3092: return "251Cf";
		case 3093: return "252Cf";
		case 3094: return "253Cf";
		case 3095: return "254Cf";
		case 3096: return "255Cf";
		case 3097: return "256Cf";
		case 3098: return "239Es";
		case 3099: return "240Es";
		case 3100: return "241Es";
		case 3101: return "242Es";
		case 3102: return "243Es";
		case 3103: return "244Es";
		case 3104: return "245Es";
		case 3105: return "246Es";
		case 3106: return "247Es";
		case 3107: return "248Es";
		case 3108: return "249Es";
		case 3109: return "250Es";
		case 3110: return "251Es";
		case 3111: return "252Es";
		case 3112: return "253Es";
		case 3113: return "254Es";
		case 3114: return "255Es";
		case 3115: return "256Es";
		case 3116: return "257Es";
		case 3117: return "258Es";
		case 3118: return "241Fm";
		case 3119: return "242Fm";
		case 3120: return "243Fm";
		case 3121: return "244Fm";
		case 3122: return "245Fm";
		case 3123: return "246Fm";
		case 3124: return "247Fm";
		case 3125: return "248Fm";
		case 3126: return "249Fm";
		case 3127: return "250Fm";
		case 3128: return "251Fm";
		case 3129: return "252Fm";
		case 3130: return "253Fm";
		case 3131: return "254Fm";
		case 3132: return "255Fm";
		case 3133: return "256Fm";
		case 3134: return "257Fm";
		case 3135: return "258Fm";
		case 3136: return "259Fm";
		case 3137: return "260Fm";
		case 3138: return "245Md";
		case 3139: return "246Md";
		case 3140: return "247Md";
		case 3141: return "248Md";
		case 3142: return "249Md";
		case 3143: return "250Md";
		case 3144: return "251Md";
		case 3145: return "252Md";
		case 3146: return "253Md";
		case 3147: return "254Md";
		case 3148: return "255Md";
		case 3149: return "256Md";
		case 3150: return "257Md";
		case 3151: return "258Md";
		case 3152: return "259Md";
		case 3153: return "260Md";
		case 3154: return "261Md";
		case 3155: return "262Md";
		case 3156: return "248No";
		case 3157: return "249No";
		case 3158: return "250No";
		case 3159: return "251No";
		case 3160: return "252No";
		case 3161: return "253No";
		case 3162: return "254No";
		case 3163: return "255No";
		case 3164: return "256No";
		case 3165: return "257No";
		case 3166: return "258No";
		case 3167: return "259No";
		case 3168: return "260No";
		case 3169: return "261No";
		case 3170: return "262No";
		case 3171: return "263No";
		case 3172: return "264No";
		case 3173: return "251Lr";
		case 3174: return "252Lr";
		case 3175: return "253Lr";
		case 3176: return "254Lr";
		case 3177: return "255Lr";
		case 3178: return "256Lr";
		case 3179: return "257Lr";
		case 3180: return "258Lr";
		case 3181: return "259Lr";
		case 3182: return "260Lr";
		case 3183: return "261Lr";
		case 3184: return "262Lr";
		case 3185: return "263Lr";
		case 3186: return "264Lr";
		case 3187: return "265Lr";
		case 3188: return "266Lr";
		case 3189: return "253Rf";
		case 3190: return "254Rf";
		case 3191: return "255Rf";
		case 3192: return "256Rf";
		case 3193: return "257Rf";
		case 3194: return "258Rf";
		case 3195: return "259Rf";
		case 3196: return "260Rf";
		case 3197: return "261Rf";
		case 3198: return "262Rf";
		case 3199: return "263Rf";
		case 3200: return "264Rf";
		case 3201: return "265Rf";
		case 3202: return "266Rf";
		case 3203: return "267Rf";
		case 3204: return "268Rf";
		case 3205: return "255Db";
		case 3206: return "256Db";
		case 3207: return "257Db";
		case 3208: return "258Db";
		case 3209: return "259Db";
		case 3210: return "260Db";
		case 3211: return "261Db";
		case 3212: return "262Db";
		case 3213: return "263Db";
		case 3214: return "264Db";
		case 3215: return "265Db";
		case 3216: return "266Db";
		case 3217: return "267Db";
		case 3218: return "268Db";
		case 3219: return "269Db";
		case 3220: return "270Db";
		case 3221: return "258Sg";
		case 3222: return "259Sg";
		case 3223: return "260Sg";
		case 3224: return "261Sg";
		case 3225: return "262Sg";
		case 3226: return "263Sg";
		case 3227: return "264Sg";
		case 3228: return "265Sg";
		case 3229: return "266Sg";
		case 3230: return "267Sg";
		case 3231: return "268Sg";
		case 3232: return "269Sg";
		case 3233: return "270Sg";
		case 3234: return "271Sg";
		case 3235: return "272Sg";
		case 3236: return "273Sg";
		case 3237: return "260Bh";
		case 3238: return "261Bh";
		case 3239: return "262Bh";
		case 3240: return "263Bh";
		case 3241: return "264Bh";
		case 3242: return "265Bh";
		case 3243: return "266Bh";
		case 3244: return "267Bh";
		case 3245: return "268Bh";
		case 3246: return "269Bh";
		case 3247: return "270Bh";
		case 3248: return "271Bh";
		case 3249: return "272Bh";
		case 3250: return "273Bh";
		case 3251: return "274Bh";
		case 3252: return "275Bh";
		case 3253: return "263Hs";
		case 3254: return "264Hs";
		case 3255: return "265Hs";
		case 3256: return "266Hs";
		case 3257: return "267Hs";
		case 3258: return "268Hs";
		case 3259: return "269Hs";
		case 3260: return "270Hs";
		case 3261: return "271Hs";
		case 3262: return "272Hs";
		case 3263: return "273Hs";
		case 3264: return "274Hs";
		case 3265: return "275Hs";
		case 3266: return "276Hs";
		case 3267: return "277Hs";
		case 3268: return "265Mt";
		case 3269: return "266Mt";
		case 3270: return "267Mt";
		case 3271: return "268Mt";
		case 3272: return "269Mt";
		case 3273: return "270Mt";
		case 3274: return "271Mt";
		case 3275: return "272Mt";
		case 3276: return "273Mt";
		case 3277: return "274Mt";
		case 3278: return "275Mt";
		case 3279: return "276Mt";
		case 3280: return "277Mt";
		case 3281: return "278Mt";
		case 3282: return "279Mt";
		case 3283: return "267Ds";
		case 3284: return "268Ds";
		case 3285: return "269Ds";
		case 3286: return "270Ds";
		case 3287: return "271Ds";
		case 3288: return "272Ds";
		case 3289: return "273Ds";
		case 3290: return "274Ds";
		case 3291: return "275Ds";
		case 3292: return "276Ds";
		case 3293: return "277Ds";
		case 3294: return "278Ds";
		case 3295: return "279Ds";
		case 3296: return "280Ds";
		case 3297: return "281Ds";
		case 3298: return "272Rg";
		case 3299: return "273Rg";
		case 3300: return "274Rg";
		case 3301: return "275Rg";
		case 3302: return "276Rg";
		case 3303: return "277Rg";
		case 3304: return "278Rg";
		case 3305: return "279Rg";
		case 3306: return "280Rg";
		case 3307: return "281Rg";
		case 3308: return "282Rg";
		case 3309: return "283Rg";
		case 3310: return "276Cn";
		case 3311: return "277Cn";
		case 3312: return "278Cn";
		case 3313: return "279Cn";
		case 3314: return "280Cn";
		case 3315: return "281Cn";
		case 3316: return "282Cn";
		case 3317: return "283Cn";
		case 3318: return "284Cn";
		case 3319: return "285Cn";
		case 3320: return "278Nh";
		case 3321: return "279Nh";
		case 3322: return "280Nh";
		case 3323: return "281Nh";
		case 3324: return "282Nh";
		case 3325: return "283Nh";
		case 3326: return "284Nh";
		case 3327: return "285Nh";
		case 3328: return "286Nh";
		case 3329: return "287Nh";
		case 3330: return "285Fl";
		case 3331: return "286Fl";
		case 3332: return "287Fl";
		case 3333: return "288Fl";
		case 3334: return "289Fl";
		case 3335: return "287Mc";
		case 3336: return "288Mc";
		case 3337: return "289Mc";
		case 3338: return "290Mc";
		case 3339: return "291Uup";
		case 3340: return "289Lv";
		case 3341: return "290Lv";
		case 3342: return "291Lv";
		case 3343: return "292Lv";
		case 3344: return "293Lv";
		case 3345: return "291Ts";
		case 3346: return "292Ts";
		case 3347: return "293Ts";
		case 3348: return "294Uus";
		case 3349: return "293Og";
		case 3350: return "294Og";
		case 3351: return "295Og";
		  default: return "?";
	}
}

/******************************************************************************

 Function nist_isotope(): return the respective enum isotope for a given atomic
 symbol s.

******************************************************************************/

isotope nist_isotope(const char s[])
{
	for (int n = 0; n < 3352; ++n)
	{
		const isotope enum_n = (isotope) n;
		if (strcmp(nist_atomic_symbol(enum_n), s) == 0) return enum_n;
	}

	PRINT_ERROR("invalid atomic symbol s = %s\n", s)
	exit(EXIT_FAILURE);
}
