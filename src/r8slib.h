int i4_max ( int i1, int i2 );
double r8_big ( );
double r8_epsilon ( );
double r8_log_2 ( double x );
double r8_max ( double x, double y );
void r8mat_add ( int m, int n, double alpha, double a[], double beta, double b[], double c[] );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
double *r8mat_copy_new ( int m, int n, double a1[] );
double *r8mat_fss_new ( int n, double a[], int nb, double b[] );
double *r8mat_identity_new ( int n );
int r8mat_is_significant ( int m, int n, double r[], double s[] );
int r8mat_minvm ( int n1, int n2, double a[], double b[], double c[] );
void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] );
double r8mat_norm_li ( int m, int n, double a[] );
void r8mat_scale ( int m, int n, double s, double a[] );
double *r8mat_zeros_new ( int m, int n );
