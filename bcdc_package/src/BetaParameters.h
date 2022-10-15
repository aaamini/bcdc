#ifndef B0B01449_7C28_4F3C_805E_47CFC4603E27
#define B0B01449_7C28_4F3C_805E_47CFC4603E27


struct BetaParameters {
    double alpha;
    double beta;

    BetaParameters() :alpha{1}, beta{1} {};
    BetaParameters(const double a, const double b) :alpha{a}, beta{b} {};

    void set(const double a, const double b)
    {
        alpha = a;
        beta = b;
    }
};


#endif /* B0B01449_7C28_4F3C_805E_47CFC4603E27 */
