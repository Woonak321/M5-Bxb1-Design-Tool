import numpy as np
import pandas as pd
from scipy.special import expit  # logistic sigmoid

class M5Model:
    """M5 Bxb1 Integration Design Tool - Core Model"""
    
    def __init__(self):
        # Mode B parameters (logit-linear)
        self.a = 3.0161  # intercept
        self.b = 1.03e-3  # slope per bp
        
        # Concentration normalization (weak dependence)
        self.K = 0.20  # μM (= 200 nM)
        self.n = 6  # Hill coefficient
        self.C_ref = 0.387  # μM (= 387 nM)
        
        # Mode A parameters
        self.C50 = 1.0  # μM
        self.nH = 1
        self.T50 = 30  # minutes
        self.T_ssAP50 = 20  # minutes
        self.t_ssAP = 45  # minutes (fixed)
        
        # DNA-level parameters
        self.k_d = 7e-4  # per kb (distance decay)
        self.gamma = 10  # GC penalty
        
        # Global calibration
        self.G = 6.632173469
        self.alpha_mech = 8.5e-4  # per bp
        self.beta_mech = 8.5e-4  # per bp
        
        # Population parameters
        self.alpha_pop = 1.0
        
        # Spacer parameters
        self.S0 = 80  # bp
        self.sigma_spacer = 40  # bp
        
        # Upper limit
        self.cap = 0.95
        
    def predict_mode_b(self, L, C=0.387, t=30):
        """
        Predict integration success rate for Mode B (production mode)
        
        Parameters:
        -----------
        L : float or array
            Insert length in kb
        C : float
            Bxb1 concentration in μM (default: 0.387)
        t : float
            Reaction time in minutes (default: 30)
            
        Returns:
        --------
        P_prod : float or array
            Predicted success probability
        """
        # Base probability (logit-linear)
        logit_P = self.a - self.b * L * 1000  # convert kb to bp
        P_base = expit(logit_P)  # sigmoid function
        
        # Concentration factor (weak dependence, reference-normalized)
        f_conc = self._concentration_factor(C)
        
        # Final probability
        P_prod = np.minimum(P_base * f_conc, self.cap)
        
        return P_prod
    
    def _concentration_factor(self, C):
        """Hill normalization relative to reference concentration"""
        numerator = C**self.n / (self.K**self.n + C**self.n)
        denominator = self.C_ref**self.n / (self.K**self.n + self.C_ref**self.n)
        return numerator / denominator
    
    def predict_mode_a(self, d, gc, L=1.0, C=0.387, t=30):
        """
        Predict integration success for Mode A (site selection)
        
        Parameters:
        -----------
        d : float
            Distance to oriC in kb
        gc : float
            GC content (0-1)
        L : float
            Insert length in kb
        C : float
            Bxb1 concentration in μM
        t : float
            Reaction time in minutes
            
        Returns:
        --------
        P_final : float
            Predicted success probability
        """
        # Distance effect
        F_distance = np.exp(-self.k_d * d)
        
        # GC penalty
        f_GC = np.exp(-self.gamma * (gc - 0.5)**2)
        
        # Protein level (phi_tot)
        phi_Bxb1 = C**self.nH / (self.C50**self.nH + C**self.nH)
        phi_time = t / (self.T50 + t)
        phi_ssAP = self.t_ssAP / (self.T_ssAP50 + self.t_ssAP)
        phi_tot = phi_Bxb1 * phi_time * phi_ssAP
        
        # Combine factors
        P_final = phi_tot * F_distance * f_GC
        
        return np.minimum(P_final, self.cap)
    
    def expected_colonies(self, L, N_CFU=30, C=0.387, alpha_pop=None):
        """
        Calculate expected positive colonies per plate
        
        Parameters:
        -----------
        L : float
            Insert length in kb
        N_CFU : int
            Colony forming units per plate
        C : float
            Bxb1 concentration in μM
        alpha_pop : float
            Population scaling factor (default: self.alpha_pop)
            
        Returns:
        --------
        expected : float
            Expected number of positive colonies
        """
        if alpha_pop is None:
            alpha_pop = self.alpha_pop
            
        P_prod = self.predict_mode_b(L, C)
        expected = N_CFU * np.minimum(alpha_pop * P_prod, 1.0)
        
        return expected
    
    def calibrate(self, observed_data, params_to_fit=['a', 'b']):
        """
        Recalibrate model with new experimental data
        
        Parameters:
        -----------
        observed_data : str or DataFrame
            Path to CSV file or pandas DataFrame with columns:
            ['insert_length_bp', 'success_rate', 'n_replicates']
        params_to_fit : list
            Parameters to optimize (default: ['a', 'b'])
            
        Returns:
        --------
        mae : float
            Mean absolute error in percentage points
        """
        from scipy.optimize import minimize
        
        if isinstance(observed_data, str):
            df = pd.read_csv(observed_data)
        else:
            df = observed_data
            
        L_obs = df['insert_length_bp'].values / 1000  # convert to kb
        P_obs = df['success_rate'].values
        
        def objective(params):
            if 'a' in params_to_fit:
                self.a = params[0]
            if 'b' in params_to_fit:
                idx = params_to_fit.index('b')
                self.b = params[idx]
                
            P_pred = self.predict_mode_b(L_obs)
            mae = np.mean(np.abs(P_pred - P_obs))
            return mae
        
        # Initial guess
        x0 = []
        if 'a' in params_to_fit:
            x0.append(self.a)
        if 'b' in params_to_fit:
            x0.append(self.b)
            
        result = minimize(objective, x0, method='Nelder-Mead')
        
        # Update parameters
        if 'a' in params_to_fit:
            self.a = result.x[0]
        if 'b' in params_to_fit:
            idx = params_to_fit.index('b')
            self.b = result.x[idx]
            
        self.mae = result.fun * 100  # convert to percentage
        
        return self.mae
    
    def save(self, filename):
        """Save model parameters to file"""
        import json
        params = {
            'a': self.a,
            'b': self.b,
            'K': self.K,
            'n': self.n,
            'C_ref': self.C_ref,
            'C50': self.C50,
            'nH': self.nH,
            'T50': self.T50,
            'k_d': self.k_d,
            'gamma': self.gamma,
            'G': self.G,
            'alpha_mech': self.alpha_mech,
            'beta_mech': self.beta_mech
        }
        with open(filename, 'w') as f:
            json.dump(params, f, indent=2)
    
    @classmethod
    def load_default(cls):
        """Load model with default parameters"""
        return cls()


if __name__ == "__main__":
    # Example usage
    model = M5Model()
    
    # Predict for 1.5 kb insert
    P = model.predict_mode_b(L=1.5)
    print(f"Predicted success rate for 1.5 kb: {P*100:.1f}%")
    
    # Calculate expected colonies
    colonies = model.expected_colonies(L=1.5, N_CFU=30)
    print(f"Expected positive colonies per plate: {colonies:.1f}")
