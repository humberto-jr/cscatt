#if !defined(GLOBALS_HEADER)
	#define GLOBALS_HEADER
	#include "c_lib.h"
	#include <omp.h>

	/******************************************************************************

	 Type real: defines a general real type number to which precision, double or
	 float, can be tuned during compilation.

	 NOTE: not used.

	******************************************************************************/

	#if defined(USE_SINGLE_PRECISION)
		typedef float real;
	#else
		typedef double real;
	#endif

	/******************************************************************************

	 Macro MAX_LINE_LENGTH(): defines a default size for strings to be tuned during
	 compilation.

	******************************************************************************/

	#if !defined(MAX_LINE_LENGTH)
		#define MAX_LINE_LENGTH 1024
	#endif

	/******************************************************************************

	 Macro INF: use an upper bound for real numbers with double precision to define
	 infinity, and likewise, -infinity.

	 NOTE: DBL_MAX from float.h is chosen, however any other reasonably large
	 number is possible.

	******************************************************************************/

	#define INF DBL_MAX

	/******************************************************************************

	 Macro M_PI(): is the value of pi as defined in the GSL library. Its usage is
	 so common that it is worth redefining it here for compatibility reasons when
	 GSL is not available.

	******************************************************************************/

	#if !defined(M_PI)
		#define M_PI 3.14159265358979323846264338328
	#endif

	/******************************************************************************

	 Macro PRINT_ERROR(): writes an error message in the C stderr using the format
	 "# file.c, function(), line n: formatted message here".

	******************************************************************************/

	#define PRINT_ERROR(format, ...)                                          \
	{                                                                         \
	  fprintf(stderr, "# %s, %s(), line %d: ", __FILE__, __func__, __LINE__); \
	  fprintf(stderr, format, ##__VA_ARGS__);                                 \
	}

	/******************************************************************************

	 Macro ASSERT(): checks the validity (true/false) of a given expression. If
	 false, it shall exit the runtime with a formatted message in the C stderr.

	******************************************************************************/

	#define ASSERT(expr)                              \
	{                                                 \
	  if (!(expr))                                    \
	  {                                               \
	    PRINT_ERROR("assertion '%s' failed\n", #expr) \
	    exit(EXIT_FAILURE);                           \
	  }                                               \
	}

	/******************************************************************************

	 Macro AS_STRING(): return a given macro as a string between double quotes.

	******************************************************************************/

	#define AS_STRING(macro) #macro

	/******************************************************************************

	 Macro PRINT_MACRO(): return the value of a given macro as a string.

	******************************************************************************/

	#define PRINT_MACRO(macro) AS_STRING(macro)

	/******************************************************************************

	 Function allocate(): allocate resources for n elements with data_size bits and
	 return a pointer pointing the first one (initialized to zero if set_zero =
	 true).

	******************************************************************************/

	static inline void *allocate(const int n,
	                             const int data_size, const bool set_zero)
	{
		void *pointer = (set_zero? calloc(n, data_size) : malloc(n*data_size));

		if (pointer == NULL)
		{
			PRINT_ERROR("unable to allocate resources for %d elements of %d bytes\n", n, data_size)
			exit(EXIT_FAILURE);
		}

		return pointer;
	}

	/******************************************************************************

	 Function as_double(): is a simple interface to a very common operation of
	 casting an integer, n, as a double precision number.

	******************************************************************************/

	inline static double as_double(const int n)
	{
		return (double) n;
	}

	/******************************************************************************

	 Function as_int(): is a simple interface to a very common operation of casting
	 an double, n, as an integer number.

	******************************************************************************/

	inline static int as_int(const double n)
	{
		return (int) n;
	}
	/******************************************************************************

	 Function wall_time(): returns the wall time in units of seconds.

	******************************************************************************/

	inline static double wall_time()
	{
		return omp_get_wtime();
	}

	/******************************************************************************

	 Function thread_id(): returns the current OpenMP thread identification.

	******************************************************************************/

	inline static int thread_id()
	{
		return omp_get_thread_num();
	}

	/******************************************************************************

	 Function max_threads(): returns the total OpenMP threads available.

	******************************************************************************/

	inline static int max_threads()
	{
		return omp_get_max_threads();
	}

	/******************************************************************************

	 Function sinc(): returns sin(x)/x for a given x.

	******************************************************************************/

	inline static double sinc(const double x)
	{
		return sin(x)/x;
	}

	/******************************************************************************

	 Function sigmoid(): return the sigmoid function at x.

	******************************************************************************/

	inline static double sigmoid(const double x)
	{
		return 1.0/(1.0 + exp(-x));
	}

	/******************************************************************************

	 Function left_trim(): trim a string s in the left side.

	******************************************************************************/

	inline static char *left_trim(char s[])
	{
		while (isspace(*s)) ++s;
		return s;
	}

	/******************************************************************************

	 Function right_trim(): trim a string s in the right side.

	******************************************************************************/

	inline static char *right_trim(char s[])
	{
		char *s_backward = s + strlen(s);

		while (isspace(*(--s_backward)));
		*(s_backward + 1) = '\0';

		return s;
	}

	/******************************************************************************

	 Function trim(): trim a string s in both sides.

	******************************************************************************/

	inline static char *trim(char s[])
	{
		return right_trim(left_trim(s));
	}

	/******************************************************************************

	 Function time_stamp(): return the current YMDHMS date as a time stamp.

	 Example: 31 May, 2001 09:45:54 AM.

	******************************************************************************/

	inline static const char *time_stamp()
	{
		time_t now = time(NULL);
		const struct tm *info = localtime(&now);

		static char stamp[50];
		strftime(stamp, 50, "%B %d, %Y %I:%M:%S %p", info);

		return stamp;
	}

	/******************************************************************************

	 Function min(): return the min between two integers a and b.

	******************************************************************************/

	inline static int min(const int a, const int b)
	{
		return (a < b? a : b);
	}

	/******************************************************************************

	 Function max(): return the max between two integers a and b.

	******************************************************************************/

	inline static int max(const int a, const int b)
	{
		return (a > b? a : b);
	}

	/******************************************************************************

	 Function factorial(): returns the factorial of n.

	 NOTE: adapted from GSL's gamma.c source code.

	******************************************************************************/

	inline static double factorial(const int n)
	{
		switch (n)
		{
			case   0: return 1.0;
			case   1: return 1.0;
			case   2: return 2.0;
			case   3: return 6.0;
			case   4: return 24.0;
			case   5: return 120.0;
			case   6: return 720.0;
			case   7: return 5040.0;
			case   8: return 40320.0;
			case   9: return 362880.0;
			case  10: return 3628800.0;
			case  11: return 39916800.0;
			case  12: return 479001600.0;
			case  13: return 6227020800.0;
			case  14: return 87178291200.0;
			case  15: return 1307674368000.0;
			case  16: return 20922789888000.0;
			case  17: return 355687428096000.0;
			case  18: return 6402373705728000.0;
			case  19: return 121645100408832000.0;
			case  20: return 2432902008176640000.0;
			case  21: return 51090942171709440000.0;
			case  22: return 1124000727777607680000.0;
			case  23: return 25852016738884976640000.0;
			case  24: return 620448401733239439360000.0;
			case  25: return 15511210043330985984000000.0;
			case  26: return 403291461126605635584000000.0;
			case  27: return 10888869450418352160768000000.0;
			case  28: return 304888344611713860501504000000.0;
			case  29: return 8841761993739701954543616000000.0;
			case  30: return 265252859812191058636308480000000.0;
			case  31: return 8222838654177922817725562880000000.0;
			case  32: return 263130836933693530167218012160000000.0;
			case  33: return 8683317618811886495518194401280000000.0;
			case  34: return 2.95232799039604140847618609644e38;
			case  35: return 1.03331479663861449296666513375e40;
			case  36: return 3.71993326789901217467999448151e41;
			case  37: return 1.37637530912263450463159795816e43;
			case  38: return 5.23022617466601111760007224100e44;
			case  39: return 2.03978820811974433586402817399e46;
			case  40: return 8.15915283247897734345611269600e47;
			case  41: return 3.34525266131638071081700620534e49;
			case  42: return 1.40500611775287989854314260624e51;
			case  43: return 6.04152630633738356373551320685e52;
			case  44: return 2.65827157478844876804362581101e54;
			case  45: return 1.19622220865480194561963161496e56;
			case  46: return 5.50262215981208894985030542880e57;
			case  47: return 2.58623241511168180642964355154e59;
			case  48: return 1.24139155925360726708622890474e61;
			case  49: return 6.08281864034267560872252163321e62;
			case  50: return 3.04140932017133780436126081661e64;
			case  51: return 1.55111875328738228022424301647e66;
			case  52: return 8.06581751709438785716606368564e67;
			case  53: return 4.27488328406002556429801375339e69;
			case  54: return 2.30843697339241380472092742683e71;
			case  55: return 1.26964033536582759259651008476e73;
			case  56: return 7.10998587804863451854045647464e74;
			case  57: return 4.05269195048772167556806019054e76;
			case  58: return 2.35056133128287857182947491052e78;
			case  59: return 1.38683118545689835737939019720e80;
			case  60: return 8.32098711274139014427634118320e81;
			case  61: return 5.07580213877224798800856812177e83;
			case  62: return 3.14699732603879375256531223550e85;
			case  63: return 1.982608315404440064116146708360e87;
			case  64: return 1.268869321858841641034333893350e89;
			case  65: return 8.247650592082470666723170306800e90;
			case  66: return 5.443449390774430640037292402480e92;
			case  67: return 3.647111091818868528824985909660e94;
			case  68: return 2.480035542436830599600990418570e96;
			case  69: return 1.711224524281413113724683388810e98;
			case  70: return 1.197857166996989179607278372170e100;
			case  71: return 8.504785885678623175211676442400e101;
			case  72: return 6.123445837688608686152407038530e103;
			case  73: return 4.470115461512684340891257138130e105;
			case  74: return 3.307885441519386412259530282210e107;
			case  75: return 2.480914081139539809194647711660e109;
			case  76: return 1.885494701666050254987932260860e111;
			case  77: return 1.451830920282858696340707840860e113;
			case  78: return 1.132428117820629783145752115870e115;
			case  79: return 8.946182130782975286851441715400e116;
			case  80: return 7.156945704626380229481153372320e118;
			case  81: return 5.797126020747367985879734231580e120;
			case  82: return 4.753643337012841748421382069890e122;
			case  83: return 3.945523969720658651189747118010e124;
			case  84: return 3.314240134565353266999387579130e126;
			case  85: return 2.817104114380550276949479442260e128;
			case  86: return 2.422709538367273238176552320340e130;
			case  87: return 2.107757298379527717213600518700e132;
			case  88: return 1.854826422573984391147968456460e134;
			case  89: return 1.650795516090846108121691926250e136;
			case  90: return 1.485715964481761497309522733620e138;
			case  91: return 1.352001527678402962551665687590e140;
			case  92: return 1.243841405464130725547532432590e142;
			case  93: return 1.156772507081641574759205162310e144;
			case  94: return 1.087366156656743080273652852570e146;
			case  95: return 1.032997848823905926259970209940e148;
			case  96: return 9.916779348709496892095714015400e149;
			case  97: return 9.619275968248211985332842594960e151;
			case  98: return 9.426890448883247745626185743100e153;
			case  99: return 9.332621544394415268169923885600e155;
			case 100: return 9.33262154439441526816992388563e157;
			case 101: return 9.42594775983835942085162312450e159;
			case 102: return 9.61446671503512660926865558700e161;
			case 103: return 9.90290071648618040754671525458e163;
			case 104: return 1.02990167451456276238485838648e166;
			case 105: return 1.08139675824029090050410130580e168;
			case 106: return 1.146280563734708354534347384148e170;
			case 107: return 1.226520203196137939351751701040e172;
			case 108: return 1.324641819451828974499891837120e174;
			case 109: return 1.443859583202493582204882102460e176;
			case 110: return 1.588245541522742940425370312710e178;
			case 111: return 1.762952551090244663872161047110e180;
			case 112: return 1.974506857221074023536820372760e182;
			case 113: return 2.231192748659813646596607021220e184;
			case 114: return 2.543559733472187557120132004190e186;
			case 115: return 2.925093693493015690688151804820e188;
			case 116: return 3.393108684451898201198256093590e190;
			case 117: return 3.96993716080872089540195962950e192;
			case 118: return 4.68452584975429065657431236281e194;
			case 119: return 5.57458576120760588132343171174e196;
			case 120: return 6.68950291344912705758811805409e198;
			case 121: return 8.09429852527344373968162284545e200;
			case 122: return 9.87504420083360136241157987140e202;
			case 123: return 1.21463043670253296757662432419e205;
			case 124: return 1.50614174151114087979501416199e207;
			case 125: return 1.88267717688892609974376770249e209;
			case 126: return 2.37217324288004688567714730514e211;
			case 127: return 3.01266001845765954480997707753e213;
			case 128: return 3.85620482362580421735677065923e215;
			case 129: return 4.97450422247728744039023415041e217;
			case 130: return 6.46685548922047367250730439554e219;
			case 131: return 8.47158069087882051098456875820e221;
			case 132: return 1.11824865119600430744996307608e224;
			case 133: return 1.48727070609068572890845089118e226;
			case 134: return 1.99294274616151887673732419418e228;
			case 135: return 2.69047270731805048359538766215e230;
			case 136: return 3.65904288195254865768972722052e232;
			case 137: return 5.01288874827499166103492629211e234;
			case 138: return 6.91778647261948849222819828311e236;
			case 139: return 9.61572319694108900419719561353e238;
			case 140: return 1.34620124757175246058760738589e241;
			case 141: return 1.89814375907617096942852641411e243;
			case 142: return 2.69536413788816277658850750804e245;
			case 143: return 3.85437071718007277052156573649e247;
			case 144: return 5.55029383273930478955105466055e249;
			case 145: return 8.04792605747199194484902925780e251;
			case 146: return 1.17499720439091082394795827164e254;
			case 147: return 1.72724589045463891120349865931e256;
			case 148: return 2.55632391787286558858117801578e258;
			case 149: return 3.80892263763056972698595524351e260;
			case 150: return 5.71338395644585459047893286526e262;
			case 151: return 8.62720977423324043162318862650e264;
			case 152: return 1.31133588568345254560672467123e267;
			case 153: return 2.00634390509568239477828874699e269;
			case 154: return 3.08976961384735088795856467036e271;
			case 155: return 4.78914290146339387633577523906e273;
			case 156: return 7.47106292628289444708380937294e275;
			case 157: return 1.17295687942641442819215807155e278;
			case 158: return 1.85327186949373479654360975305e280;
			case 159: return 2.94670227249503832650433950735e282;
			case 160: return 4.71472363599206132240694321176e284;
			case 161: return 7.59070505394721872907517857094e286;
			case 162: return 1.22969421873944943411017892849e289;
			case 163: return 2.00440157654530257759959165344e291;
			case 164: return 3.28721858553429622726333031164e293;
			case 165: return 5.42391066613158877498449501421e295;
			case 166: return 9.00369170577843736647426172359e297;
			case 167: return 1.50361651486499904020120170784e300;
			case 168: return 2.52607574497319838753801886917e302;
			case 169: return 4.26906800900470527493925188890e304;
			case 170: return 7.25741561530799896739672821113e306;

			/* NOTE: larger n-values may cause overflow. */
			default: return as_double(n)*factorial(n - 1);
		}
	}
#endif
