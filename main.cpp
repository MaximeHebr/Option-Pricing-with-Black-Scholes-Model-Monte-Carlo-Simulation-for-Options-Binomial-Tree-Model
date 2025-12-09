#include "Options.h"
#include "Market.h"
#include "Pricers.h"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;

void printHeader(const string& title) {
    cout << "\n" << string(70, '=') << "\n";
    cout << title << "\n";
    cout << string(70, '=') << "\n";
}

template<typename Func>
double measureTime(Func&& func) {
    auto start = chrono::high_resolution_clock::now();
    func();
    auto end = chrono::high_resolution_clock::now();
    return chrono::duration<double, milli>(end - start).count();
}

void example1_BasicComparison() {
    printHeader("Exemple 1: Comparaison des trois méthodes");
    
    auto option = make_shared<VanillaOption>(100.0, 1.0, OptionType::CALL);
    auto market = make_shared<MarketData>(100.0, 0.05, 0.20, 0.0);
    
    cout << "\nOption: Call européen, Strike=100, Maturité=1 an\n";
    cout << "Marché: S=100, r=5%, σ=20%\n\n";
    
    // Black-Scholes
    BlackScholesPricer bs(option, market);
    double timeBS = measureTime([&]() { bs.calculatePrice(); });
    auto result = bs.calculateFullResult();
    
    cout << "BLACK-SCHOLES:\n";
    cout << result.toString() << "\n";
    cout << "Temps: " << fixed << setprecision(2) << timeBS << " ms\n\n";
    
    // Binomial
    BinomialTreePricer binomial(option, market, 1000);
    double timeBT = measureTime([&]() { binomial.calculatePrice(); });
    double priceB = binomial.calculatePrice();
    
    cout << "ARBRE BINOMIAL (1000 steps):\n";
    cout << "Prix: " << setprecision(4) << priceB << "\n";
    cout << "Erreur: " << abs(priceB - result.price) << "\n";
    cout << "Temps: " << setprecision(2) << timeBT << " ms\n\n";
    
    // Monte Carlo
    MonteCarloPricer mc(option, market, 100000);
    double timeMC = measureTime([&]() { mc.calculatePrice(); });
    auto [priceMC, errorMC] = mc.calculatePriceWithError();
    
    cout << "MONTE CARLO (100k simulations):\n";
    cout << "Prix: " << setprecision(4) << priceMC << " ± " << errorMC << "\n";
    cout << "Erreur: " << abs(priceMC - result.price) << "\n";
    cout << "Temps: " << setprecision(2) << timeMC << " ms\n";
}

void example2_AmericanOption() {
    printHeader("Exemple 2: Option américaine vs européenne");
    
    auto american = make_shared<VanillaOption>(110.0, 1.0, OptionType::PUT, ExerciseStyle::AMERICAN);
    auto european = make_shared<VanillaOption>(110.0, 1.0, OptionType::PUT);
    auto market = make_shared<MarketData>(90.0, 0.05, 0.20, 0.0);
    
    cout << "\nPut profondément dans la monnaie (S=90, K=110)\n\n";
    
    BlackScholesPricer bsEuro(european, market);
    double euroPrice = bsEuro.calculatePrice();
    cout << "Put européen (BS): " << fixed << setprecision(4) << euroPrice << "\n";
    
    BinomialTreePricer btAmer(american, market, 1000);
    double amerPrice = btAmer.calculatePrice();
    cout << "Put américain (Binomial): " << amerPrice << "\n";
    
    cout << "\nPrime d'exercice anticipé: " << amerPrice - euroPrice << "\n";
}

void example3_AsianOption() {
    printHeader("Exemple 3: Option asiatique");
    
    auto asian = make_shared<AsianOption>(100.0, 1.0, OptionType::CALL);
    auto vanilla = make_shared<VanillaOption>(100.0, 1.0, OptionType::CALL);
    auto market = make_shared<MarketData>(100.0, 0.05, 0.25, 0.0);
    
    cout << "\nComparaison Call vanille vs Call asiatique\n\n";
    
    BlackScholesPricer bsVanilla(vanilla, market);
    double vanillaPrice = bsVanilla.calculatePrice();
    cout << "Call vanille: " << fixed << setprecision(4) << vanillaPrice << "\n";
    
    MonteCarloPricer mcAsian(asian, market, 50000, 252);
    auto [asianPrice, asianError] = mcAsian.calculatePriceWithError();
    cout << "Call asiatique: " << asianPrice << " ± " << asianError << "\n";
    
    cout << "\nL'option asiatique est moins chère (volatilité réduite)\n";
    cout << "Différence: " << vanillaPrice - asianPrice << "\n";
}

void runAllExamples() {
    example1_BasicComparison();
    cout << "\nAppuyez sur Entrée pour continuer...";
    cin.ignore();
    cin.get();
    
    example2_AmericanOption();
    cout << "\nAppuyez sur Entrée pour continuer...";
    cin.ignore();
    
    example3_AsianOption();
}

int main() {
    cout << fixed << setprecision(4);
    
    while (true) {
        cout << "\n" << string(70, '=') << "\n";
        cout << "FRAMEWORK DE PRICING D'OPTIONS\n";
        cout << string(70, '=') << "\n";
        cout << "\n1. Comparaison des trois méthodes\n";
        cout << "2. Option américaine\n";
        cout << "3. Option asiatique\n";
        cout << "4. Tous les exemples\n";
        cout << "0. Quitter\n";
        cout << "\nChoix: ";
        
        int choice;
        cin >> choice;
        
        if (cin.fail()) {
            cin.clear();
            cin.ignore(10000, '\n');
            continue;
        }
        
        try {
            switch (choice) {
                case 0: return 0;
                case 1: example1_BasicComparison(); break;
                case 2: example2_AmericanOption(); break;
                case 3: example3_AsianOption(); break;
                case 4: runAllExamples(); break;
                default: cout << "\nChoix invalide\n";
            }
        } catch (const exception& e) {
            cout << "\nErreur: " << e.what() << "\n";
        }
        
        if (choice != 4) {
            cout << "\nAppuyez sur Entrée pour continuer...";
            cin.ignore(10000, '\n');
            cin.get();
        }
    }
    
    return 0;
}
