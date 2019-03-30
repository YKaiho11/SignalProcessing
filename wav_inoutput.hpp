#pragma once

typedef struct
{
    int fs; //サンプリング周波数
    int bits; //量子化bit数
    int L; //データ長
} WAV_PRM;

double *audio_read(WAV_PRM *prm, char *filename)
{
    //変数宣言
    FILE *fp;
    int n;
    double *data;
    char header_ID[4];
    long header_size;
    char header_type[4];
    char fmt_ID[4];
    long fmt_size;
    short fmt_format;
    short fmt_channel;
    long fmt_samples_per_sec;
    long fmt_bytes_per_sec;
    short fmt_block_size;
    short fmt_bits_per_sample;
    char data_ID[4];
    long data_size;
    short data_data;
    
    //wavファイルオープン
    fp = fopen(filename, "rb");
    
    //wavデータ読み込み
    fread(header_ID, 1, 4, fp);
    fread(&header_size, 4, 1, fp);
    fread(header_type, 1, 4, fp);
    fread(fmt_ID, 1, 4, fp);
    fread(&fmt_size, 4, 1, fp);
    fread(&fmt_format, 2, 1, fp);
    fread(&fmt_channel, 2, 1, fp);
    fread(&fmt_samples_per_sec, 4, 1, fp);
    fread(&fmt_bytes_per_sec, 4, 1, fp);
    fread(&fmt_block_size, 2, 1, fp);
    fread(&fmt_bits_per_sample, 2, 1, fp);
    char c[2];
    for(int i=0;i<100;i++){
        fread(c,1,2,fp);
        printf("%s\n",c);
        if(!strcmp(c,"da")){
            fread(c,1,2,fp);
            printf("%s\n",c);
            if(!strcmp(c,"ta")) break;
        }
    }
    fread(&data_size, 4, 1, fp);
    
    //パラメータ代入
    prm->fs = fmt_samples_per_sec;
    prm->bits = fmt_bits_per_sample;
    prm->L = (int)data_size / 2;
    
    //音声データ代入
    data = (double*)calloc(prm->L,sizeof(double));
    for (n = 0; n < prm->L; n++) {
        fread(&data_data, 2, 1, fp);
        data[n] = (double)data_data / 32768.0;
    }
    
    fclose(fp);
    return data;
}


void audio_write(double *data, WAV_PRM *prm, char *filename)
{
    //変数宣言
    FILE *fp;
    int n;
    char header_ID[4];
    long header_size;
    char header_type[4];
    char fmt_ID[4];
    long fmt_size;
    short fmt_format;
    short fmt_channel;
    long fmt_samples_per_sec;
    long fmt_bytes_per_sec;
    short fmt_block_size;
    short fmt_bits_per_sample;
    char data_ID[4];
    long data_size;
    short data_data;
    
    //ファイルオープン
    fp = fopen(filename, "w");
    
    //ヘッダー書き込み
    header_ID[0] = 'R';
    header_ID[1] = 'I';
    header_ID[2] = 'F';
    header_ID[3] = 'F';
    header_size = 36 + prm->L * 2;
    header_type[0] = 'W';
    header_type[1] = 'A';
    header_type[2] = 'V';
    header_type[3] = 'E';
    fwrite(header_ID, 1, 4, fp);
    fwrite(&header_size, 4, 1, fp);
    fwrite(header_type, 1, 4, fp);
    
    //フォーマット書き込み
    fmt_ID[0] = 'f';
    fmt_ID[1] = 'm';
    fmt_ID[2] = 't';
    fmt_ID[3] = ' ';
    fmt_size = 16;
    fmt_format = 1;
    fmt_channel = 1;
    fmt_samples_per_sec = prm->fs;
    fmt_bytes_per_sec = prm->fs * prm->bits / 8;
    fmt_block_size = prm->bits / 8;
    fmt_bits_per_sample = prm->bits;
    fwrite(fmt_ID, 1, 4, fp);
    fwrite(&fmt_size, 4, 1, fp);
    fwrite(&fmt_format, 2, 1, fp);
    fwrite(&fmt_channel, 2, 1, fp);
    fwrite(&fmt_samples_per_sec, 4, 1, fp);
    fwrite(&fmt_bytes_per_sec, 4, 1, fp);
    fwrite(&fmt_block_size, 2, 1, fp);
    fwrite(&fmt_bits_per_sample, 2, 1, fp);
    
    //データ書き込み
    data_ID[0] = 'd';
    data_ID[1] = 'a';
    data_ID[2] = 't';
    data_ID[3] = 'a';
    data_size = prm->L * 2;
    fwrite(data_ID, 1, 4, fp);
    fwrite(&data_size, 4, 1, fp);
    
    //音声データ書き込み
    //fp = fopen(filename, "w");
    if(fp==NULL) printf("failed to open file!");
    for (n = 0; n < prm->L; n++) {
        //リミッター
        if (data[n] > 1) {
            data_data = 32767;
        } else if (data[n] < -1) {
            data_data = -32767;
        } else {
            data_data = (short)(data[n] * 32767.0);
        }
        fwrite(&data_data, 2, 1, fp);
    }
    
    fclose(fp);
}
