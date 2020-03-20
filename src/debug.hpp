/* Include this hpp file in any file where you need to output
 * cerr or cout. The advantage of using this file is that if you
 * just comment out the "define DEBUG" line, then all the places
 * where you used "op" to output cout's and "err" to output cerr's,
 * will become inactive. ie: They won't output anything. Keeps the
 * code clean and efficient.

 * Usage: In your code, instead of typing std::cout<<"Hi";
 * Just type cout("Hi");
 * Same way, instead of std::cerr<<"Error"; just type cerr("Error");
 */

#ifndef DEBUGHELPER_HPP
#define DEBUGHELPER_HPP

#define DEBUG_MODE

#define USEFUL_CODE_TEMPORARILY_COMMENTED_OUT_DONT_DELETE_IT
#define CODE_FOR_TESTING_FEEL_FREE_TO_DELETE_IT

#ifdef DEBUG_MODE
    //Colour code by Vaughan Schmidt. http://www.codebuilder.me/2014/01/color-terminal-text-in-c/
    #define RESET   "\033[37m"      // previously \033[0m
    #define BLACK   "\033[30m"      // Black
    #define RED     "\033[31m"      // Red
    #define GREEN   "\033[32m"      // Green
    #define YELLOW  "\033[33m"      // Yellow
    #define BLUE    "\033[34m"      // Blue
    #define MAGENTA "\033[35m"      // Magenta
    #define CYAN    "\033[36m"      // Cyan
    #define WHITE   "\033[37m"      // White
    #define BOLDBLACK   "\033[1m\033[30m"      // Bold Black
    #define BOLDRED     "\033[1m\033[31m"      // Bold Red
    #define BOLDGREEN   "\033[1m\033[32m"      // Bold Green
    #define BOLDYELLOW  "\033[1m\033[33m"      // Bold Yellow
    #define BOLDBLUE    "\033[1m\033[34m"      // Bold Blue
    #define BOLDMAGENTA "\033[1m\033[35m"      // Bold Magenta
    #define BOLDCYAN    "\033[1m\033[36m"      // Bold Cyan
    #define BOLDWHITE   "\033[1m\033[37m"      // Bold White
    #define CLEAR "\033[2J"  // clear screen escape code
    #define coutBlack(x) (std::cout << BLACK << (x) << RESET)
    #define coutBoldBlack(x) (std::cout << BOLDBLACK << (x) << RESET)
    #define coutRed(x) (std::cout << RED << (x) << RESET)
    #define coutBoldRed(x) (std::cout << BOLDRED << (x) << RESET)
    #define coutGreen(x) (std::cout << GREEN << (x) << RESET)
    #define coutBoldGreen(x) (std::cout << BOLDGREEN << (x) << RESET)
    #define coutYellow(x) (std::cout << YELLOW << (x) << RESET)
    #define coutBoldYellow(x) (std::cout << BOLDYELLOW << (x) << RESET)
    #define coutBlue(x) (std::cout << BLUE << (x) << RESET)
    #define coutBoldBlue(x) (std::cout << BOLDBLUE << (x) << RESET)
    #define coutMagenta(x) (std::cout << MAGENTA << (x) << RESET)
    #define coutBoldMagenta(x) (std::cout << BOLDMAGENTA << (x) << RESET)
    #define coutCyan(x) (std::cout << CYAN << (x) << RESET)
    #define coutBoldCyan(x) (std::cout << BOLDCYAN << (x) << RESET)
    #define coutWhite(x) (std::cout << WHITE << (x) << RESET)
    #define coutBoldWhite(x) (std::cout << BOLDWHITE << (x) << RESET)
    #define cerr(x) (std::cerr << (x))
    #define cout(x) (std::cout << (x))
    //... etc
    #define Q() (std::exit(0))

#else
    #define cerr(x)
    #define cout(x)
    //... etc


#endif

#endif // DEBUGHELPER_HPP
