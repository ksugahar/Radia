// Read-only CUDA API observation in an owned process; never call a CUDA API.
// Frida module observer contract: https://frida.re/docs/javascript-api/
const hooked=new Set();
const context=/^(cuCtxCreate(?:_v[0-9]+)?|cuDevicePrimaryCtxRetain|cuCtxSetCurrent|cuCtxPushCurrent(?:_v[0-9]+)?|cudaSetDevice|cudaMalloc.*|cudaFree)$/;
const observed=/^(cuInit|cuGetProcAddress(?:_v[0-9]+)?|cuGetExportTable|cuDeviceGetCount|cudaGetDeviceCount|cudaGetDriverEntryPoint.*)$/;
let sequence=0;
function emit(value){value.sequence=++sequence;send(value);}
function hook(address,name,origin){
    const key=address.toString();if(hooked.has(key))return;hooked.add(key);
    try{
        Interceptor.attach(address,{
            onEnter(args){
                this.symbol=null;this.output=null;
                if(/^(cuGetProcAddress|cudaGetDriverEntryPoint)/.test(name)){
                    try{this.symbol=args[0].readUtf8String();this.output=args[1];}
                    catch(error){emit({kind:'coverage_error',name:name,error:String(error)});}
                }
                emit({kind:context.test(name)?'context_api':'observed_api',name:name,origin:origin,symbol:this.symbol,thread:Process.getCurrentThreadId()});
            },
            onLeave(result){
                if(this.symbol && this.output && result.toInt32()===0 && (context.test(this.symbol)||observed.test(this.symbol))){
                    try{const target=this.output.readPointer();if(!target.isNull())hook(target,this.symbol,'dynamic-resolution');
                        else emit({kind:'coverage_error',name:this.symbol,error:'null successful entrypoint'});}
                    catch(error){emit({kind:'coverage_error',name:this.symbol,error:String(error)});}
                }
            }
        });
        emit({kind:'hook',name:name,origin:origin});
    }catch(error){emit({kind:'coverage_error',name:name,error:String(error)});}
}
const observer=Process.attachModuleObserver({
    onAdded(module){
        if(/nvcuda|cudart|radia_cuda_csr_pcg/i.test(module.name)){
            emit({kind:'module',name:module.name,path:module.path});
            if(/nvcuda|cudart/i.test(module.name)){
                try{for(const item of module.enumerateExports()){
                    if(item.type==='function' && (context.test(item.name)||observed.test(item.name)))hook(item.address,item.name,module.name);
                }}catch(error){emit({kind:'coverage_error',name:module.name,error:String(error)});}
            }
        }
    }
});
emit({kind:'ready',pid:Process.id});
