#include "treeutils.hpp"

int main(){
    string s = "((((((((((((((((((TAEGU:0.042517,GEOFO:0.036656)1.000000:0.066334,CORBR:0.033967)1.000000:0.071515,MANVI:0.081597)1.000000:0.023207,ACACH:0.106500)1.000000:0.052243,(MELUN:0.068663,NESNO:0.045636)1.000000:0.051793)1.000000:0.003149,FALPE:0.097938)1.000000:0.004553,CARCR:0.067141)1.000000:0.004408,(((((PICPU:0.199928,MERNU:0.140702)0.992950:0.004658,BUCRH:0.132149)1.000000:0.006483,APAVI:0.134081)0.999794:0.004106,LEPDI:0.088891)1.000000:0.005531,COLST:0.161849)0.981887:0.003333)0.961321:0.001695,(((HALAL:0.006139,HALLE:0.006194)1.000000:0.046530,CATAU:0.035196)1.000000:0.005933,TYTAL:0.087578)0.998797:0.001288)1.000000:0.008523,(CHAVO:0.079614,BALRE:0.075611)0.998466:0.002263)0.536554:0.000186,(((((((EGRGA:0.070778,PELCR:0.052009)1.000000:0.004303,NIPNI:0.051311)1.000000:0.005794,PHACA:0.083001)1.000000:0.007294,((PYGAD:0.012740,APTFO:0.009195)1.000000:0.038119,FULGL:0.041360)1.000000:0.002796)1.000000:0.006869,GAVST:0.054400)1.000000:0.004917,(PHALE:0.071653,EURHE:0.131300)1.000000:0.007244)0.647633:0.001247,OPHHO:0.097300)0.424196:0.000190)0.438111:0.000357,((CHAPE:0.103963,CALAN:0.136451)1.000000:0.061501,CAPCA:0.086462)1.000000:0.008538)0.526540:0.001528,((TAUER:0.091219,CHLUN:0.091803)0.998352:0.001732,CUCCA:0.142809)0.999038:0.003587)0.969343:0.001459,(((PODCR:0.076110,PHORU:0.033127)1.000000:0.042716,COLLI:0.121842)0.503534:0.000723,(MESUN:0.110408,PTEGU:0.087200)1.000000:0.006877)0.995149:0.001703)1.000000:0.028403,(TINMA:0.178284,STRCA:0.132000)1.000000:0.074699)1.000000:0.029143,ANAPL:0.144782)1.000000:0.008008,MELGA:0.064782),GALGA);";
    Tree t(s);
    cerr << t.newick() << endl;
    BinaryTree bt(t);
    cerr << bt.newick() << endl;
    vector<BinaryTree> bts = BinaryTree::text2trees(s + s + s);
    cerr << bts.size() << endl;
}