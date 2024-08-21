#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
#include <RcppGSL.h>
#include <gsl/gsl_sf_psi.h> 

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
double T_density(double x,
                 double Mu,
                 double Sigma,
                 double Nu){
  return (Rf_gammafn((Nu + 1) / 2.0) / (Rf_gammafn(Nu / 2.0)* sqrt(M_PI * Nu * pow(Sigma, 2))))*pow((1+pow((x-Mu),2)/(Nu*pow(Sigma,2))),-(Nu+1)/2.0);
}

// [[Rcpp::export]]
double Mix_T_density_x(double x,
                     vec Pi, 
                     vec Mu,
                     vec Sigma,
                     vec Nu){
  vec dens_of_x_for_each_t(Pi.n_elem);
  for (int i =0; i< Pi.n_elem;++i){
    dens_of_x_for_each_t[i]=Pi[i]*T_density(x,Mu[i],Sigma[i],Nu[i]);
    }
  double densdens_of_X=sum(dens_of_x_for_each_t);
  return densdens_of_X;
}



// [[Rcpp::export]]
vec Mix_T_density(vec X,
                     vec Pi, 
                     vec Mu,
                     vec Sigma,
                     vec Nu){
  vec Dens_for_each_x(X.n_elem);
  
  for (int j=0; j<X.n_elem;++j){
    Dens_for_each_x[j] = Mix_T_density_x(X[j],Pi,Mu,Sigma,Nu);
  }
  return (Dens_for_each_x);
}


// [[Rcpp::export]]
double likelyhood_t_mix(vec X,
                  vec Pi,
                  vec mu,
                  vec Sigma,
                  vec Nu){
  return prod(Mix_T_density(X,Pi,mu,Sigma,Nu));
}

// [[Rcpp::export]]
double loglikelyhood_t_mix(vec X,
                        vec Pi,
                        vec mu,
                        vec Sigma,
                        vec Nu){
  vec densities = Mix_T_density(X, Pi, mu, Sigma, Nu);
  double log_likelihood = sum(log(densities));
  return log_likelihood;
}


// [[Rcpp::export]]
mat Weights_of_X(vec X,
                 vec Pi,
                 vec mu,
                 vec Sigma,
                 vec Nu){
  //Row's is a given outcome Pi is the probiallty that it come from the j density
  mat Weigths(X.n_elem,Pi.n_elem);
  for (int i =0; i < X.n_elem; ++i){
    vec like_x_i(Pi.n_elem); // make vector
    for (int j=0;j<Pi.n_elem;++j){
      like_x_i[j]=Pi[j]*T_density(X[i],mu[j],Sigma[j],Nu[j]);
    }
    for (int k=0;k<Pi.n_elem;++k){
    Weigths(i,k)=like_x_i[k]/sum(like_x_i);
    }
  }
  return Weigths;
}



// [[Rcpp::export]]
vec get_Pi(mat Weigths){
  if (Weigths.n_rows == 0) {
    Rcpp::stop("Input matrix has no rows");
  }
  return sum(Weigths, 0).t() / Weigths.n_rows;
}



// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
mat gradient(vec X,
             vec Pi,
             vec mu,
             vec Sigma,
             vec Nu,
             string normalise ="no_normaliztion"){
  // col mu ,Sigma, Nu
  mat Gradiets(Pi.n_elem,3);
  // find weights
  mat w_x = Weights_of_X(X,Pi,mu,Sigma,Nu);
  
  
  // mu
  for (int i =0; i < mu.n_elem; ++i){
    double mu_i =0;
    for (int j =0; j < X.n_elem; ++j){
      // mising weight
      mu_i=mu_i+w_x(j,i)*(X[j]-mu[i])/(Nu[i]*pow(Sigma[i],2)+pow(X[j]-mu[i],2));
    }
    mu_i = mu_i*(Nu[i]+1);
    mu_i =-mu_i; 
    Gradiets(i,0)=mu_i;
  // sigma
    double sigma_i=0;
    for (int j =0; j < X.n_elem; ++j){
      sigma_i= sigma_i+w_x(j,i)*(-1+(Nu[i]+1)*pow(X[j]-mu[i],2)/(Sigma[i]*Nu[i]+pow(X[j]-mu[i],2)));
      }
    sigma_i=sigma_i/Sigma[i];
    sigma_i=-sigma_i;
    Gradiets(i,1)=sigma_i;
    // nu
    // declare varibels for less coput
    double Gamma_v_plus_one= tgamma((Nu[i]+1)/2);
    double psi_v_plus_one= gsl_sf_psi((Nu[i]+1)/2);
    double Gamma_v = tgamma(Nu[i]/2);
    double psi_v = gsl_sf_psi(Nu[i]/2);
    double nu_i=0;
    for (int j =0;j < X.n_elem;++j){
      nu_i=nu_i-0.5*w_x(j,i)+((Nu[i]+1)/2)*((pow(X[j]-mu[i],2))/(pow(Sigma[i]*Nu[i],2)))*w_x(j,i)+((psi_v_plus_one)/(2*Gamma_v_plus_one))*w_x(j,i)-((psi_v)/(2*Gamma_v))*w_x(j,i);
      }
    nu_i=-nu_i;
    Gradiets(i,2)=nu_i;
  }
    // Gradient normalsation
  if (normalise =="second_norm"){
    Gradiets=Gradiets/norm(Gradiets,2); // second option based on dataX.n_elem;
    }
  if (normalise =="N_element"){
    Gradiets=Gradiets/X.n_elem;
    }
  if (normalise =="Element_wise_Normalsation"){
    Gradiets=Gradiets.col(0)/norm(Gradiets.col(0),2);
    Gradiets=Gradiets.col(1)/norm(Gradiets.col(1),2);
    Gradiets=Gradiets.col(2)/norm(Gradiets.col(2),2);
  }
  return Gradiets;
}

// [[Rcpp::export]]
vec Maximum_of_two_vector(vec v1,vec v2){
  vec v_out(v1.n_elem);
  for(int i =0;i<v1.n_elem;++i){
    v_out[i]= max(v1[i],v2[i]);
  }
  return v_out;
}


//using overload to define to function based on the case where in

mat Gradient_clipping_vec(mat gradient,vec C_levels=vec(3,1.0)){
  for (int i =0;i<gradient.n_cols;++i){
    if(norm(gradient.col(i),2)>C_levels[i]){
      gradient.col(i)=C_levels[i]*gradient.col(i)/norm(gradient.col(i),2);
      }
    }
  return(gradient);
}
mat Gradient_clipping(mat gradient,double C_level=1){
  if(norm(gradient,2)>C_level){
    gradient=C_level*gradient/norm(gradient,2);
  }
  return(gradient);
}
// Chosen not to export graident clibbing to R since i cant export to function with the same name.

/*
I this section i have made multiple version of gradient decent.
It could be combinede in to one function. I have keept it as threes versions so, since in the final EM one is gona be called as the standart.
So insted of investing af lot time to writing it, this makes the flow faster.
*/
// [[Rcpp::export]]
mat gradient_decent_no_check(vec X,
                                   vec Pi,
                                   vec mu,
                                   vec Sigma,
                                   vec Nu,
                                   double alpha=0.02,
                                   double Max_iter=1000,
                                   double norm_gradient = 0.5,
                                   string normalise_method ="no_normaliztion"//this should be changed
                               )
  {
  mat parameters(Pi.n_elem,3); 
  // Set parameters value in matrix form (esay to update)
  parameters.col(0)=mu;
  parameters.col(1)=Sigma;
  parameters.col(2)=Nu;
  mat gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),normalise_method); 
  
  // Sigma_lower_bound_vec
  /*
  This is smart since it the if one overshoot the domain of the Gamma function is not definede 
  */
  double parameter_min= 0.01;
  vec Sigma_lower_bound_vec(Sigma.n_elem,fill::ones);
  Sigma_lower_bound_vec=parameter_min*Sigma_lower_bound_vec;
  vec Nu_lower_bound_vec(Nu.n_elem,fill::ones);
  Nu_lower_bound_vec=parameter_min*Nu_lower_bound_vec;
  
  // loop throug gradient decent
  for (int i=0;i<Max_iter;++i){
    if (norm(vectorise(gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),"no_normaliztion")),2)<norm_gradient){
      return parameters;
    }
    parameters=parameters-alpha*gradient_in_point;
    // If undershoot set to minimum Note is is importante that the check is done before the gradient is updated, otherwise the gradient will not be able to be evaluated
    parameters.col(1)=Maximum_of_two_vector(parameters.col(1),Sigma_lower_bound_vec);
    parameters.col(2)=Maximum_of_two_vector(parameters.col(2),Nu_lower_bound_vec);
    // update parametes
    gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),normalise_method);
    
  }
  return parameters;
}


// [[Rcpp::export]]
mat gradient_decent_with_check(vec X,
                             vec Pi,
                             vec mu,
                             vec Sigma,
                             vec Nu,
                             double alpha=0.02,
                             double Max_iter=1000,
                             double norm_gradient = 0.5,
                             string normalise_method ="no_normaliztion",//this should be changed
                             double max_it_in_check= 10
                            )
{
  mat parameters(Pi.n_elem,3); 
  // Set parameters value in matrix form (esay to update)
  parameters.col(0)=mu;
  parameters.col(1)=Sigma;
  parameters.col(2)=Nu;
  mat gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),normalise_method); 
  // make teperay variabel
  mat tem_parameters=parameters;
  // Sigma_lower_bound_ve
  double parameter_min= 0.01;
  vec Sigma_lower_bound_vec(Sigma.n_elem,fill::ones);
  Sigma_lower_bound_vec=parameter_min*Sigma_lower_bound_vec;
  vec Nu_lower_bound_vec(Nu.n_elem,fill::ones);
  Nu_lower_bound_vec=parameter_min*Nu_lower_bound_vec;
  // loop throug gradient decent
  for (int i=0;i<Max_iter;++i){
    if (norm(vectorise(gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),"no_normaliztion")),2)<norm_gradient){
      return parameters;
    }
    for(int j=0;j<max_it_in_check;j++){
      if(likelyhood_t_mix(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2))<likelyhood_t_mix(X,Pi,tem_parameters.col(0),tem_parameters.col(1),tem_parameters.col(2))){
        break;
      }
      tem_parameters=parameters-alpha*gradient_in_point;
      gradient_in_point=alpha*gradient_in_point;
    }
    parameters=tem_parameters;
    // If undershoot set to minimum Note is is importante that the check is done before the gradient is updated, otherwise the gradient will not be able to be evaluated
    parameters.col(1)=Maximum_of_two_vector(parameters.col(1),Sigma_lower_bound_vec);
    parameters.col(2)=Maximum_of_two_vector(parameters.col(2),Nu_lower_bound_vec);
    // update parametes
    gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),normalise_method);
    
  }
  return parameters;
}

// need to be made
// [[Rcpp::export]]
mat gradient_decent_T_mix_model(vec X,
                               vec Pi,
                               vec mu,
                               vec Sigma,
                               vec Nu,
                               double alpha=0.02,
                               double Max_iter=1000,
                               double norm_gradient = 0.5,
                               string type_decent ="Normal",
                               string normalise_method ="no_normaliztion",//this should be changed
                               double max_it_in_check= 10
                               
)
{
  mat parameters(Pi.n_elem,3); 
  // Set parameters value in matrix form (esay to update)
  parameters.col(0)=mu;
  parameters.col(1)=Sigma;
  parameters.col(2)=Nu;
  mat gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),normalise_method); 
  // make teperay variabel
  mat tem_parameters=parameters;
  // Sigma_lower_bound_ve
  double parameter_min= 0.01;
  vec Sigma_lower_bound_vec(Sigma.n_elem,fill::ones);
  Sigma_lower_bound_vec=parameter_min*Sigma_lower_bound_vec;
  vec Nu_lower_bound_vec(Nu.n_elem,fill::ones);
  Nu_lower_bound_vec=parameter_min*Nu_lower_bound_vec;
  // loop throug gradient decent
  for (int i=0;i<Max_iter;++i){
    if (norm(vectorise(gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),"no_normaliztion")),2)<norm_gradient){
      return parameters;
    }
    for(int j=0;j<max_it_in_check;j++){
      if(likelyhood_t_mix(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2))<likelyhood_t_mix(X,Pi,tem_parameters.col(0),tem_parameters.col(1),tem_parameters.col(2))){
        break;
      }
      tem_parameters=parameters-alpha*gradient_in_point;
      gradient_in_point=alpha*gradient_in_point;
    }
    parameters=tem_parameters;
    // If undershoot set to minimum Note is is importante that the check is done before the gradient is updated, otherwise the gradient will not be able to be evaluated
    parameters.col(1)=Maximum_of_two_vector(parameters.col(1),Sigma_lower_bound_vec);
    parameters.col(2)=Maximum_of_two_vector(parameters.col(2),Nu_lower_bound_vec);
    // update parametes
    gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),normalise_method);
    
  }
  return parameters;
}

// Gradient decent whit clipping
//
// [[Rcpp::export]]
mat gradient_decent_with_clipping(vec X,
                                vec Pi,
                                vec mu,
                                vec Sigma,
                                vec Nu,
                                vec Clipping_vector,
                                double alpha=0.02,
                                int Max_iter=1000,
                                double norm_gradient = 0.5,
                                int max_it_in_check= 10
)
{
  // make mat for parameters
  mat parameters(Pi.n_elem,3); 
  
  // Set parameters value in matrix form (esay to update)
  parameters.col(0)=mu;
  parameters.col(1)=Sigma;
  parameters.col(2)=Nu;
  mat gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),"no_normaliztion"); 
  
  // make teperay variabel
  mat tem_parameters=parameters;
  
  // Make lover bound for parameter
  double parameter_min= 0.01;
  vec Sigma_lower_bound_vec(Sigma.n_elem,fill::ones);
  Sigma_lower_bound_vec=parameter_min*Sigma_lower_bound_vec;
  vec Nu_lower_bound_vec(Nu.n_elem,fill::ones);
  Nu_lower_bound_vec=parameter_min*Nu_lower_bound_vec;
  
  // loop throug gradient decent
  for (int i=0;i<Max_iter;++i){
    // make check for gradient size 
    if (norm(vectorise(gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),"no_normaliztion")),2)<norm_gradient){
      return parameters;
    }
    //graident clipping prcedorure
    //mat Gradient_clipping_vec(mat gradient,vec C_levels=vec(3,1.0))
    gradient_in_point=Gradient_clipping_vec(gradient_in_point,Clipping_vector);
    
    //check if gradient result in small likelyhood.
    for(int j=0;j<max_it_in_check;j++){
      if(likelyhood_t_mix(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2))<likelyhood_t_mix(X,Pi,tem_parameters.col(0),tem_parameters.col(1),tem_parameters.col(2))){
        break;
      }
      tem_parameters=parameters-alpha*gradient_in_point;
      gradient_in_point=alpha*gradient_in_point; // this contiues aplies alpha to the graiden 
    }
    parameters=tem_parameters;
    // If undershoot set to minimum Note is is importante that the check is done before the gradient is updated, otherwise the gradient will not be able to be evaluated
    parameters.col(1)=Maximum_of_two_vector(parameters.col(1),Sigma_lower_bound_vec);
    parameters.col(2)=Maximum_of_two_vector(parameters.col(2),Nu_lower_bound_vec);
    // update parametes
    gradient_in_point = gradient(X,Pi,parameters.col(0),parameters.col(1),parameters.col(2),"no_normaliztion");
    
  }
  return parameters;
}






// [[Rcpp::export]]
mat optimize_using_gradient_decent(vec X,
                                   vec Pi,
                                   vec mu,
                                   vec Sigma,
                                   vec Nu,
                                   double alpha=0.02,
                                   int Max_iter=1000,
                                   double norm_gradient = 0.5,
                                   string normalise_method ="no_normaliztion",//this should be changed
                                   bool Check_gradient_decrease = true,
                                   int max_it_in_check= 100){
  if(Check_gradient_decrease==false){
    return(gradient_decent_no_check(X,Pi,mu,Sigma,Nu,alpha,Max_iter,norm_gradient,normalise_method));
  }
  else{
    return(gradient_decent_with_check(X,Pi,mu,Sigma,Nu,alpha,Max_iter,norm_gradient,normalise_method,max_it_in_check));
  }
  
}


// [[Rcpp::export]]
mat EM_t_distribution(vec X,
                      vec Pi,
                      vec mu,
                      vec Sigma,
                      vec Nu,
                      vec Clipping_vector,
                      int Max_iter = 10000,
                      double alpha=0.02,
                      int max_iter_gradient=100,
                      double norm_gradient = 0.001
){
  mat parameters(Pi.n_elem,4); 
  // Set parameters value in matrix form (esay to update)
  
  parameters.col(0)=mu;
  parameters.col(1)=Sigma;
  parameters.col(2)=Nu;
  
  mat w_x=Weights_of_X(X,Pi,mu,Sigma,Nu);
  parameters.col(3)=get_Pi(w_x);
  for (int i=0;i<Max_iter;++i){
    w_x=Weights_of_X(X,parameters.col(3),parameters.col(0),parameters.col(1),parameters.col(2));
    parameters.col(3)=get_Pi(w_x);
    parameters.cols(0,2)=gradient_decent_with_clipping(X,parameters.col(3),parameters.col(0),parameters.col(1),parameters.col(2),Clipping_vector,alpha,Max_iter,norm_gradient,max_iter_gradient);
  }
  
  return parameters;
}


