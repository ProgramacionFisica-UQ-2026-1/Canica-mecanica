#importando las librerias 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from abc import ABC, abstractmethod
#
#
class Sistema(ABC):
    """Clase base abstracta para sistemas dinámicos."""
    
    @abstractmethod
    def derivada(self, t, estado):
        """Retorna el vector de derivadas dado el tiempo y el estado."""
        pass
    
    @abstractmethod
    def simular(self, t_span, estado_inicial, t_eval=None):
        """Simula el sistema y retorna tiempos y solución."""
        pass
#
#
class Lorenz(Sistema):
    """Modelo del atractor de Lorenz."""
    
    def __init__(self, sigma=10.0, rho=28.0, beta=8/3):
        self.sigma = sigma
        self.rho = rho
        self.beta = beta
    
    def derivada(self, t, estado):
        x, y, z = estado
        dxdt = self.sigma * (y - x)
        dydt = x * (self.rho - z) - y
        dzdt = x * y - self.beta * z
        return [dxdt, dydt, dzdt]
    
    def simular(self, t_span, estado_inicial, t_eval=None):
        sol = solve_ivp(self.derivada, t_span, estado_inicial,
                        t_eval=t_eval, method='RK45')
        return sol.t, sol.y  # sol.y tiene forma (3, N)

#
#
class Analizador:
    """Genera datos sintéticos, ajusta parámetros y visualiza el atractor de Lorenz."""
    
    def __init__(self, modelo_real, estado_inicial, t_max, n_puntos, ruido_std=0.5):
        self.modelo = modelo_real          # instancia de Lorenz con parámetros "reales"
        self.estado0 = np.array(estado_inicial)
        self.t_max = t_max
        self.n_puntos = n_puntos
        self.ruido_std = ruido_std
        
        # Tiempos de medición equiespaciados (donde simularemos "medir")
        self.t_med = np.linspace(0, t_max, n_puntos)
        
        # Simular con el modelo real y añadir ruido solo a la variable x (la que mediremos)
        _, y_real = modelo_real.simular((0, t_max), self.estado0, t_eval=self.t_med)
        self.x_real = y_real[0]  # primera fila = x(t)
        ruido = np.random.normal(0, ruido_std, size=self.x_real.shape)
        self.x_med = self.x_real + ruido   # datos sintéticos de x
        
        # Para guardar resultados del ajuste
        self.param_ajustados = None

  #
  #
      def modelo_para_ajuste(self, t, sigma, rho, beta):
        # Creamos una instancia temporal con estos parámetros
        temp_model = Lorenz(sigma, rho, beta)
        _, y_sim = temp_model.simular((0, self.t_max), self.estado0, t_eval=t)
        # Devolvemos solo la componente x (primera fila)
        return y_sim[0]
#
#funcion para los parametros
    def ajustar_parametros(self):
        # Valores iniciales para sigma, rho, beta (distintos de los reales)
        p0 = [5.0, 20.0, 2.0]
        # curve_fit minimiza la diferencia entre x_med y modelo_para_ajuste
        popt, pcov = curve_fit(self.modelo_para_ajuste, self.t_med, self.x_med,
                               p0=p0, maxfev=10000)
        self.param_ajustados = popt
        self.incertidumbre = np.sqrt(np.diag(pcov))
        return popt, pcov

#Graficas 
    def graficar_resultados(self):
        if self.param_ajustados is None:
            self.ajustar_parametros()
        
        # Crear un modelo con los parámetros ajustados para simular la curva suave
        sigma_aj, rho_aj, beta_aj = self.param_ajustados
        modelo_ajustado = Lorenz(sigma_aj, rho_aj, beta_aj)
        
        # Tiempos densos para curvas suaves
        t_denso = np.linspace(0, self.t_max, 500)
        _, y_real_denso = self.modelo.simular((0, self.t_max), self.estado0, t_eval=t_denso)
        _, y_ajust_denso = modelo_ajustado.simular((0, self.t_max), self.estado0, t_eval=t_denso)
        
        # Figura 1: comparación de la serie temporal de x
        plt.figure(figsize=(12,4))
        plt.plot(t_denso, y_real_denso[0], 'b-', label='x (modelo real)')
        plt.scatter(self.t_med, self.x_med, color='red', s=10, label='Datos sintéticos')
        plt.plot(t_denso, y_ajust_denso[0], 'g--', label='x (modelo ajustado)')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.legend()
        plt.title('Componente x: real vs datos vs ajustado')
        plt.grid(True)
        plt.show()
        
        # Figura 2: atractor 3D
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(y_real_denso[0], y_real_denso[1], y_real_denso[2], 'b-', linewidth=0.8)
        ax.set_title('Atractor de Lorenz (parámetros reales)')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()
        
        # Figura 3: comparación en el plano XZ (proyección)
        plt.figure()
        plt.plot(y_real_denso[0], y_real_denso[2], 'b-', label='Modelo real')
        plt.plot(y_ajust_denso[0], y_ajust_denso[2], 'g--', label='Modelo ajustado')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.legend()
        plt.title('Proyección XZ')
        plt.grid(True)
        plt.show()

  #analisis de errrores     def analisis_error(self):
        if self.param_ajustados is None:
            self.ajustar_parametros()
        x_pred = self.modelo_para_ajuste(self.t_med, *self.param_ajustados)
        errores = self.x_med - x_pred
        mse = np.mean(errores**2)
        mae = np.mean(np.abs(errores))
        print("----- ANÁLISIS DE ERROR -----")
        print(f"Error cuadrático medio: {mse:.4f}")
        print(f"Error absoluto medio: {mae:.4f}")
        print(f"Ruido inyectado (std): {self.ruido_std}")
        return errores

  if __name__ == "__main__":
    # Parámetros "reales" (los que intentaremos recuperar)
    sigma_real = 10.0
    rho_real = 28.0
    beta_real = 8/3
    estado0 = [1.0, 1.0, 1.0]   # condición inicial típica
    t_max = 30                    # tiempo de simulación (en las unidades de Lorenz)
    n_puntos = 200                # número de puntos "medidos"
    ruido = 0.8                   # desviación estándar del ruido
    
    # Crear modelo real
    lorenz_real = Lorenz(sigma_real, rho_real, beta_real)
    
    # Crear analizador
    analizador = Analizador(lorenz_real, estado0, t_max, n_puntos, ruido)
    
    # Ajustar parámetros
    popt, pcov = analizador.ajustar_parametros()
    print("Parámetros ajustados (sigma, rho, beta):", popt)
    print("Incertidumbre (1 sigma):", analizador.incertidumbre)
    
    # Análisis de error
    analizador.analisis_error()
    
    # Graficar todo
    analizador.graficar_resultados()
