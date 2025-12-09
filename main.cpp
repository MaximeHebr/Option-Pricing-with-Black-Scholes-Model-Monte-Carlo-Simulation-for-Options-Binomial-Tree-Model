#include "Options.h"
#include "Market.h"
#include "Pricers.h"
#include <iostream>
#include <iomanip>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#endif

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
    
    double spot, strike, maturity, rate, vol, dividend;
    char optionTypeChoice;
    
    cout << "\n=== PARAMÈTRES DE L'OPTION ===\n";
    cout << "Prix spot (ex: 100): ";
    cin >> spot;
    cout << "Strike (ex: 100): ";
    cin >> strike;
    cout << "Maturité en années (ex: 1): ";
    cin >> maturity;
    cout << "Type (C pour Call, P pour Put): ";
    cin >> optionTypeChoice;
    
    cout << "\n=== PARAMÈTRES DE MARCHÉ ===\n";
    cout << "Taux sans risque en % (ex: 5): ";
    cin >> rate;
    rate /= 100.0;  // Convert to decimal
    cout << "Volatilité en % (ex: 20): ";
    cin >> vol;
    vol /= 100.0;  // Convert to decimal
    cout << "Dividende en % (ex: 0): ";
    cin >> dividend;
    dividend /= 100.0;  // Convert to decimal
    
    OptionType type = (optionTypeChoice == 'C' || optionTypeChoice == 'c') ? OptionType::CALL : OptionType::PUT;
    
    auto option = make_shared<VanillaOption>(strike, maturity, type);
    auto market = make_shared<MarketData>(spot, rate, vol, dividend);
    
    cout << "\nOption: " << (type == OptionType::CALL ? "Call" : "Put") 
         << " européen, Strike=" << strike << ", Maturité=" << maturity << " an(s)\n";
    cout << "Marché: S=" << spot << ", r=" << (rate*100) << "%, σ=" << (vol*100) << "%\n\n";
    
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
    
    double spot, strike, maturity, rate, vol;
    
    cout << "\n=== PARAMÈTRES DE L'OPTION PUT ===\n";
    cout << "Prix spot (ex: 90): ";
    cin >> spot;
    cout << "Strike (ex: 110): ";
    cin >> strike;
    cout << "Maturité en années (ex: 1): ";
    cin >> maturity;
    
    cout << "\n=== PARAMÈTRES DE MARCHÉ ===\n";
    cout << "Taux sans risque en % (ex: 5): ";
    cin >> rate;
    rate /= 100.0;
    cout << "Volatilité en % (ex: 20): ";
    cin >> vol;
    vol /= 100.0;
    
    auto american = make_shared<VanillaOption>(strike, maturity, OptionType::PUT, ExerciseStyle::AMERICAN);
    auto european = make_shared<VanillaOption>(strike, maturity, OptionType::PUT);
    auto market = make_shared<MarketData>(spot, rate, vol, 0.0);
    
    cout << "\nPut (S=" << spot << ", K=" << strike << ")\n\n";
    
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
    
    double spot, strike, maturity, rate, vol;
    char optionTypeChoice;
    
    cout << "\n=== PARAMÈTRES DE L'OPTION ===\n";
    cout << "Prix spot (ex: 100): ";
    cin >> spot;
    cout << "Strike (ex: 100): ";
    cin >> strike;
    cout << "Maturité en années (ex: 1): ";
    cin >> maturity;
    cout << "Type (C pour Call, P pour Put): ";
    cin >> optionTypeChoice;
    
    cout << "\n=== PARAMÈTRES DE MARCHÉ ===\n";
    cout << "Taux sans risque en % (ex: 5): ";
    cin >> rate;
    rate /= 100.0;
    cout << "Volatilité en % (ex: 25): ";
    cin >> vol;
    vol /= 100.0;
    
    OptionType type = (optionTypeChoice == 'C' || optionTypeChoice == 'c') ? OptionType::CALL : OptionType::PUT;
    
    auto asian = make_shared<AsianOption>(strike, maturity, type);
    auto vanilla = make_shared<VanillaOption>(strike, maturity, type);
    auto market = make_shared<MarketData>(spot, rate, vol, 0.0);
    
    string typeStr = (type == OptionType::CALL) ? "Call" : "Put";
    cout << "\nComparaison " << typeStr << " vanille vs " << typeStr << " asiatique\n\n";
    
    BlackScholesPricer bsVanilla(vanilla, market);
    double vanillaPrice = bsVanilla.calculatePrice();
    cout << typeStr << " vanille: " << fixed << setprecision(4) << vanillaPrice << "\n";
    
    MonteCarloPricer mcAsian(asian, market, 50000, 252);
    auto [asianPrice, asianError] = mcAsian.calculatePriceWithError();
    cout << typeStr << " asiatique: " << asianPrice << " ± " << asianError << "\n";
    
    cout << "\nL'option asiatique est généralement moins chère (volatilité réduite)\n";
    cout << "Différence: " << abs(vanillaPrice - asianPrice) << "\n";
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
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif
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
                case 1: 
                    cin.ignore(10000, '\n');
                    example1_BasicComparison(); 
                    break;
                case 2: 
                    cin.ignore(10000, '\n');
                    example2_AmericanOption(); 
                    break;
                case 3: 
                    cin.ignore(10000, '\n');
                    example3_AsianOption(); 
                    break;
                case 4: 
                    cin.ignore(10000, '\n');
                    runAllExamples(); 
                    break;
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
