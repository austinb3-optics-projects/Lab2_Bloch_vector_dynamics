# Notes on alternate scripts


- A large portion of the final script was completed by the author: Austin Burrington
- These alternate scripts were generated by AI mostly to format images
- They did also account for an inconsistency in the time evolution operator implementation that the author originally used and what was more accurate.
    - That being said, the author is still unclear how his original approach differs from the new approach.




## Time evolution operator

- A potentail reason why the new implementation was more accurate was because it utilized matrix multiplication for each product corresponding to a new time step.
- The formula for the time evolution operator used in this script is: $$\hat{\boldsymbol{U}}(t,t_{0}) \approx \prod^{N-1}_{n=0}\exp\left\{-\frac{i}{\hbar}\delta t \hat{H}(t_{0} + n\delta t)\right\}$$

- The formula (not sure if it will render) essentially allows us to account for time evolution over short time steps if the hamiltonian does not commute with itself at different times.

- The original implementation used `U_temp(:,:,i) = expm(-1i./hbar.*dt.*H);` and then did `U(:,:,i) = prod(U_temp,3)` as the mechanism to perform the product.
- the new method used `U_temp = expm(-1i/hbar*dt*H);` and `U(:,:,i) = U_temp* U(:,:,i-1)`
