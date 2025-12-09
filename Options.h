#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>

enum class OptionType { CALL, PUT };
enum class ExerciseStyle { EUROPEAN, AMERICAN };

// Base class for all options
class Option {
protected:
    double strike_;
    double maturity_;
    OptionType type_;
    ExerciseStyle style_;
    
public:
    Option(double strike, double maturity, OptionType type, ExerciseStyle style = ExerciseStyle::EUROPEAN)
        : strike_(strike), maturity_(maturity), type_(type), style_(style) {
        if (strike <= 0 || maturity <= 0) {
            throw std::invalid_argument("Strike and maturity must be positive"); // on vérifie que le strike et la maturité sont positifs pour que ce soit cohérent avec le contexte financier.
        }
    }
    
    virtual ~Option() = default;
    virtual double payoff(double spotPrice) const = 0;
    // Getters pour éviter l'accès direct aux membres privés
    double getStrike() const { return strike_; }
    double getMaturity() const { return maturity_; }
    OptionType getType() const { return type_; }
    ExerciseStyle getStyle() const { return style_; }

    bool isCall() const { return type_ == OptionType::CALL; }
    bool isPut() const { return type_ == OptionType::PUT; }
    bool isEuropean() const { return style_ == ExerciseStyle::EUROPEAN; }
    bool isAmerican() const { return style_ == ExerciseStyle::AMERICAN; }
};

// Vanilla European/American options
class VanillaOption : public Option {
public:
    VanillaOption(double strike, double maturity, OptionType type, 
                  ExerciseStyle style = ExerciseStyle::EUROPEAN)
        : Option(strike, maturity, type, style) {}
    
    double payoff(double spotPrice) const override {
        if (isCall()) {
            return std::max(spotPrice - strike_, 0.0);
        } else {
            return std::max(strike_ - spotPrice, 0.0);
        }
    }
};

// Asian option (path-dependent)
// Une option asiatique dont le payoff dépend du prix moyen du sous-jacent sur une période donnée
class AsianOption : public Option {
private:
    std::vector<double> pricePath_;
    
public:
    AsianOption(double strike, double maturity, OptionType type)
        : Option(strike, maturity, type, ExerciseStyle::EUROPEAN) {}
    
    void setPricePath(const std::vector<double>& path) { pricePath_ = path; }
    
    double calculateAverage() const {
        if (pricePath_.empty()) return 0.0;
        return std::accumulate(pricePath_.begin(), pricePath_.end(), 0.0) / pricePath_.size();
    }
    
    double payoff(double spotPrice) const override {
        double avgPrice = pricePath_.empty() ? spotPrice : calculateAverage();
        
        if (isCall()) {
            return std::max(avgPrice - strike_, 0.0);
        } else {
            return std::max(strike_ - avgPrice, 0.0);
        }
    }
};

#endif
