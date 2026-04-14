#include <lal/LALSimIMR.h>

#include <lal/LALDict.h>

#include <lal/LALValue.h>

#include <LALSimIMRPhenomX.h>
#include <LALSimIMRPhenomX_internals.h>

#include <stdio.h>
#include <stdlib.h>


void print_htilde22_properties(const COMPLEX16FrequencySeries *htilde) {
    if (!htilde) {
        printf("htilde22 is NULL\n");
        return;
    }

    // 输出基本属性
    printf("========= htilde22 Properties =========\n");
    printf("Name: %s\n", htilde->name ? htilde->name : "[unnamed]");
    printf("Start Frequency (f0): %.9f Hz\n", htilde->f0);
    printf("Frequency Spacing (deltaF): %.9e Hz\n", htilde->deltaF);

    // 输出时间戳
    char gps_str[64];
    XLALGPSToStr(gps_str, &htilde->epoch);
    printf("Epoch: %s\n", gps_str);

    // // 输出单位信息
    // char unit_str[64];
    // XLALUnitToString(unit_str);
    // printf("Units: %s\n", unit_str);

    // 输出数据信息
    if (htilde->data) {
        const UINT4 show_num = (htilde->data->length > 10) ? 10 : htilde->data->length;
        printf("Total Data Points: %u\n", htilde->data->length);
        printf("First %u Complex Samples:\n", show_num);
        
        for (UINT4 i = 0; i < show_num; ++i) {
            const REAL8 freq = htilde->f0 + i * htilde->deltaF;
            const COMPLEX16 val = htilde->data->data[i];
            printf("  [%04u] Freq: %9.3f Hz \t Real: %+.6e \t Imag: %+.6ei\n",
                   i, freq, creal(val), cimag(val));
        }
    } else {
        printf("Data: NULL\n");
    }
    printf("=======================================\n");
}

int main() {
    COMPLEX16FrequencySeries *htilde22 = NULL;

    const REAL8Sequence *freqs;
    IMRPhenomXWaveformStruct *pWF;

    // REAL8* m1_SI_ptr = malloc(sizeof(REAL8));
    // *m1_SI_ptr = 30.0;
    // const REAL8* m1_SI = m1_SI_ptr;
    // free(m1_SI_ptr);

    const REAL8 m1_SI = 30.0*LAL_MSUN_SI;
    const REAL8 m2_SI = 40.0*LAL_MSUN_SI;
    const REAL8 chi1L = 0.0;     const REAL8 chi2L = 0.0;
    const REAL8 f_min = 15.0;    const REAL8 f_max = 2048.0; REAL8 deltaF = 1.0/16.0;
    const REAL8 phiRef = 0.0;    const REAL8 fRef_in = 0;
    const REAL8 distance = 100.0 * 1e6 * LAL_PC_SI;

    LALDict* lalParams;
    lalParams = XLALCreateDict();

    // XLALDictInsertREAL8Value(lalParams, "deltaF", 1/64.0);

    // XLALPrintErrorMessage;
    // XLALSimIMRPhenomXPHMOneMode;
    // IMRPhenomXASGenerateFD;

    XLALSimIMRPhenomXASGenerateFD(&htilde22, m1_SI, m2_SI, chi1L, chi2L, distance, f_min, f_max, deltaF, phiRef, fRef_in, lalParams);

    printf("123");

    // 在main函数中调用
    print_htilde22_properties(htilde22);

    FILE *fp = fopen("htilde22.txt", "w");
    // if (!fp || !htilde22 || !htilde22->data) {
    //     XLAL_ERROR("文件或数据指针无效");
    // }

    printf("123");

    const REAL8 f0 = htilde22->f0;         // 起始频率 (Hz)
    // const REAL8 deltaF = htilde22->deltaF; // 频率间隔
    const UINT4 length = htilde22->data->length;

    for (UINT4 i = 0; i < length; ++i) {
        const REAL8 freq = f0 + i * deltaF;
        const COMPLEX16 val = htilde22->data->data[i];
        fprintf(fp, "%.9f\t%.9e\t%.9e\n", 
                freq, 
                creal(val),  // 实部 (real part)
                cimag(val)); // 虚部 (imaginary part)
    }
    printf("123");
    fclose(fp);


    if (htilde22) XLALDestroyCOMPLEX16FrequencySeries(htilde22);
    if (lalParams) XLALDestroyDict(lalParams);
    return 0;
}

