from scipy.optimize import fsolve
import math

def calculate_uncertainty(equation, theta_initial, variable_names, best_estimates, uncertainties):
  """
  Calculates the uncertainty in the solution to an equation using the numerical method.

  Args:
    equation: The equation to be solved, as a callable function.
    variable_names: A list of the variable names in the equation.
    best_estimates: A list of the best estimates for each variable.
    uncertainties: A list of the uncertainties for each variable.

  Returns:
    The estimated uncertainty in the solution.
  """

  # Initialize an empty list to store deviations
  deviations = []

  # Loop through each variable
  for i in range(len(variable_names)):
      # Calculate the modified value by adding the uncertainty to the best estimate
      modified_value = best_estimates[i] + uncertainties[i]

      # Create a new set of arguments for fsolve using the modified value
      modified_args = best_estimates[:]
      modified_args[i] = modified_value

      # Use fsolve to find the solution with the modified arguments
      theta_modified = fsolve(equation, theta_initial, args=tuple(modified_args))

      # Calculate the absolute difference between theta_modified and theta_initial and append it to the deviations list
      deviation = abs(theta_modified[0] - theta_initial)
      deviations.append(deviation)

  # Return the maximum value in the deviations list as the estimated uncertainty
  return max(deviations)

def equation_to_solve(theta, *args):
  """
  The equation to be solved, rearranged to equal zero.

  Args:
    theta: The ramp angle in radians.
    args: A list of other variables: Ft, g, mc, mb, mu_k.

  Returns:
    The result of the equation.
  """
  Ft, g, mc, mb, mu_k = args
  return Ft - g * (mc + mb) * math.sin(theta) - mu_k * mb * g * math.cos(theta)

def calculate_u_theta_k(hyp, h, u_hyp, u_h):
  """
  Calculates the uncertainty of theta_k (kinetic friction angle)

  Args:
    hyp: Hypotenuse of triangle with degree theta_k
    h: Opposite leg of triangle with degree theta_k (height)
    u_hyp: Uncertainty in hypotenuse
    u_h: Uncertainty in height
  """
  return (1 / (hyp * math.sqrt(1 - (h**2 / hyp**2)))) * math.sqrt(u_h**2 + u_hyp**2 * (h**2 / hyp**2))

if __name__ == "__main__":
  # Constants
  g = 9.8  # Acceleration due to gravity (m/s^2)
  theta = 3.6733 * (math.pi / 180) # The calculated angle (in radians) of the ramp using best estimates
  hypotenuse = 1.589 # Length of the ramp or the hypotenuse of a triangle with angle theta_k (m) 
  height = 0.438 # Elevation of the ramp of height of the opposite side of theta_k (m)

  # Best estimates for variables
  mc = 0.31112  # Mass of cart (kg)
  mb = 0.03195  # Mass of box (kg)
  theta_k = 16 * math.pi / 180  # Kinetic friction angle (radians)
  mu_k = math.tan(theta_k)  # Coefficient of kinetic friction
  Ft = 0.3051  # Thrust force (N)

  # Uncertainties
  u_g = 0.2 # Uncertainty of gravity assumed to be 0.2 (m/s^2)
  u_mc = 0.000003  # Uncertainty in mass of cart (kg)
  u_mb = 0.000003  # Uncertainty in mass of box (kg)
  u_hyp = 0.0004 # Uncertainty in hypotenuse of triangle with angle theta_k (m)
  u_h = 0.0004 # Uncertainty in hypotenuse of triangle with angle theta_k (m)
  u_theta_k = calculate_u_theta_k(hypotenuse, height, u_hyp, u_h)  # Uncertainty in kinetic friction angle (radians)
  u_mu_k = (1 / math.cos(theta_k))**2 * u_theta_k  # Uncertainty in coefficient of kinetic friction
  u_Ft = 0.01531  # Uncertainty in thrust force (N)

  # Variable names and their corresponding best estimates and uncertainties
  variable_names = ["Ft", "g", "mc", "mb", "mu_k"]
  best_estimates = [Ft, g, mc, mb, mu_k]  # Initial guess for theta, then best estimates
  uncertainties = [u_Ft, u_g, u_mc, u_mb, u_mu_k]  # Include uncertainty for g

  # Calculate the uncertainty in theta
  uncertainty_theta = calculate_uncertainty(
      lambda theta, Ft, g, mc, mb, mu_k: equation_to_solve(theta, Ft, g, mc, mb, mu_k),
      theta,
      variable_names,
      best_estimates,
      uncertainties
  )

  # Print the result
  print("Uncertainty in theta:", uncertainty_theta)
  print("Uncertainty in theta (degrees):", uncertainty_theta * 180 / math.pi)