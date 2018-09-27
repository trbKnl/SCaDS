
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List sparseSCAcpp( arma::mat X, int Q, double RIDGE, arma::vec LASSO, 
                       arma::mat fixW, int maxItrOuterloop, int nStarts,
                       bool print, double tol){ 
                       

    // Initialize objects
    
    Rcpp::List ret;

    arma::mat W0;
    arma::mat W;
    arma::mat P;
    arma::mat XTX = X.t() * X;
    arma::vec XTXdiag = XTX.diag();
    arma::vec lossFunctionValue;
    arma::vec loss;
    arma::vec sumPRESS;
    int J = X.n_cols;
    int I = X.n_rows;
    double wold;
    double CP;
    int sign_CP;
    double wols;
    double wnew;
    double minLoss; 
    double minLossGlobal = arma::datum::inf; 
 
     arma::mat U2;
     arma::vec D2; 
     arma::mat V2;
 
    
    //Initialize iterators
    int z;
    int n;
    int q;
    int j;

    //Initialize W
    arma::mat U;
    arma::vec D;
    arma::mat V;
    arma::svd(U, D, V, X);
    W0 = V.cols(0, (Q - 1));
    W0.elem( find(fixW == 0 ) ).zeros();
    W = W0;


    arma::mat test;
    // Multi-Start Cycle
    for( z = 0; z < nStarts; z++){

        lossFunctionValue = arma::vec( maxItrOuterloop + 1 );
        lossFunctionValue.fill( arma::datum::inf );

        for(n = 0; n < maxItrOuterloop; n++){

             // Update P
             arma::svd(U2, D2, V2, (XTX * W));
             P = U2.cols(0, (Q-1)) * V2.t();

             // Update W: Coordinate Descent
             for( q = 0; q < Q; q++ ){
                
                sumPRESS =  X * ( P.col(q) - W.col(q) );
                
                for( j = 0; j < J; j++ ){
                    if( fixW(j, q) != 0 ){
                        wold = W(j, q);
                        CP  = (1 / double(I)) * (dot((X.col(j)),
                                    sumPRESS.t()) + (wold * XTXdiag(j)));  
                        sign_CP = (CP > 0) ? 1 : ((CP < 0) ? -1 : 0);
                        wols = double(sign_CP) * (std::abs(CP) - LASSO(q));
                        wnew = wols / ( RIDGE + XTXdiag(j) * (1/double(I)));

                        if( std::abs(CP) < LASSO(q) ) {

                            wnew = 0;
                            sumPRESS += X.col(j) * wold;

                        } else {

                            sumPRESS += (wold - wnew) * X.col(j);
                        }

                        W(j, q) = wnew;
                    }

                }

            }
            // End Coordinate Descent
            // Update loss function 
            lossFunctionValue( n ) = ( 1 / double(I * 2) ) * accu( pow( X - ( X * W * P.t() ) , 2 ) ) + ( (1 / double(2)) * RIDGE * accu( pow( W, 2) ) ) + dot(LASSO, sum( abs(W), 0 )); 
            if( print ){
                Rcpp::Rcout << lossFunctionValue( n ) << "\n";
            } 
            
            // Evaluate condition to stop inner loop 
            if( n > 0 && lossFunctionValue( n - 1 ) - lossFunctionValue( n ) <= tol ){
                Rcpp::Rcout << "converged" << "\n";
                break; 
            }
        
        } 

        loss = lossFunctionValue( find(lossFunctionValue != arma::datum::inf) );
        minLoss = loss( loss.size() - 1 ); 
 
        if( minLoss < minLossGlobal ){
            
            minLossGlobal = minLoss; 
            
            ret["W"] = W; 
            ret["P"] = P; 
            ret["loss"] = minLoss; 
        }   
 
    }
    
    return ret;
}




