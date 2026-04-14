// lalsimulation_hook.js
console.log("[*] Starting LALSimulation hook script...");

// 确保在目标库加载后执行Hook
function hook_lalsimulation() {
    // 指定目标库的名称，根据你的环境调整
    var target_library = "lalsimulation";
    
    // 等待目标库加载
    Module.ensureInitialized(target_library).then(function() {
        console.log("[+] " + target_library + " library loaded");
        
        // Hook IMRPhenomXGetPhaseCoefficients
        var IMRPhenomXGetPhaseCoefficients_ptr = Module.findExportByName(target_library, "IMRPhenomXGetPhaseCoefficients");
        if (IMRPhenomXGetPhaseCoefficients_ptr) {
            console.log("[+] Found IMRPhenomXGetPhaseCoefficients at: " + IMRPhenomXGetPhaseCoefficients_ptr);
            Interceptor.attach(IMRPhenomXGetPhaseCoefficients_ptr, {
                onEnter: function(args) {
                    console.log("\n=== IMRPhenomXGetPhaseCoefficients Called ===");
                    // 打印参数
                    console.log("pWF: " + args[0]);
                    console.log("pPhase: " + args[1]);
                    // 在这里添加你的自定义逻辑
                },
                onLeave: function(retval) {
                    console.log("IMRPhenomXGetPhaseCoefficients return: " + retval);
                    // 如果需要，可以修改返回值，例如：retval.replace(0);
                }
            });
        } else {
            console.log("[-] IMRPhenomXGetPhaseCoefficients not found");
        }
        
        // Hook IMRPhenomX_Inspiral_Phase_22_Ansatz
        var IMRPhenomX_Inspiral_Phase_22_Ansatz_ptr = Module.findExportByName(target_library, "IMRPhenomX_Inspiral_Phase_22_Ansatz");
        if (IMRPhenomX_Inspiral_Phase_22_Ansatz_ptr) {
            console.log("[+] Found IMRPhenomX_Inspiral_Phase_22_Ansatz at: " + IMRPhenomX_Inspiral_Phase_22_Ansatz_ptr);
            Interceptor.attach(IMRPhenomX_Inspiral_Phase_22_Ansatz_ptr, {
                onEnter: function(args) {
                    console.log("\n=== IMRPhenomX_Inspiral_Phase_22_Ansatz Called ===");
                    console.log("Mf: " + args[0].toDouble());
                    console.log("powers_of_Mf: " + args[1]);
                    console.log("pPhase: " + args[2]);
                },
                onLeave: function(retval) {
                    console.log("IMRPhenomX_Inspiral_Phase_22_Ansatz return: " + retval.toDouble());
                    // 示例：修改返回值
                    // var newValue = 1.5;
                    // retval.replace(ptr(newValue));
                }
            });
        }
        
        // Hook IMRPhenomXASGenerateFD
        var IMRPhenomXASGenerateFD_ptr = Module.findExportByName(target_library, "IMRPhenomXASGenerateFD");
        if (IMRPhenomXASGenerateFD_ptr) {
            console.log("[+] Found IMRPhenomXASGenerateFD at: " + IMRPhenomXASGenerateFD_ptr);
            Interceptor.attach(IMRPhenomXASGenerateFD_ptr, {
                onEnter: function(args) {
                    console.log("\n=== IMRPhenomXASGenerateFD Called ===");
                    console.log("htilde22: " + args[0]);
                    console.log("freqs_In: " + args[1]);
                    console.log("pWF: " + args[2]);
                    console.log("lalParams: " + args[3]);
                },
                onLeave: function(retval) {
                    console.log("IMRPhenomXASGenerateFD return: " + retval);
                }
            });
        }

        // 类似的，你可以为其他函数添加Hook
        // IMRPhenomX_Inspiral_Phase_22_AnsatzInt
        // IMRPhenomX_Phase_22_ConnectionCoefficients
        // IMRPhenomX_Phase_22
        // IMRPhenomX_dPhase_22
        // IMRPhenomXPGenerateFD
        
    }).catch(function(err) {
        console.log("[-] Error waiting for library: " + err.message);
    });
}

// 使用dlopen Hook确保在库加载时执行
function hook_dlopen() {
    var android_dlopen_ext = Module.findExportByName(null, "android_dlopen_ext");
    if (android_dlopen_ext) {
        Interceptor.attach(android_dlopen_ext, {
            onEnter: function(args) {
                var so_name = args[0].readCString();
                if (so_name && so_name.indexOf("lalsimulation") >= 0) {
                    this.call_hook = true;
                    console.log("[*] lalsimulation library being loaded: " + so_name);
                }
            },
            onLeave: function(retval) {
                if (this.call_hook) {
                    setTimeout(hook_lalsimulation, 100);
                }
            }
        });
    }
}

// 主函数
function main() {
    // 尝试立即Hook，如果库已加载
    hook_lalsimulation();
    
    // 同时设置dlopen Hook以捕获后续加载
    hook_dlopen();
}

// 启动脚本
setImmediate(main);