#define __FUNC__
#define INFO(msg) PRINT('("INFO : [",A,"](",A,") : ",A)'),__FILE__,__FUNC__,msg
#define WARNING(msg) PRINT('("WARNING : [",A,"](",A,") : ",A)'),__FILE__,__FUNC__,msg
#define ERROR(msg) PRINT('("ERROR : [",A,"](",A,") : ",A)'),__FILE__,__FUNC__,msg
