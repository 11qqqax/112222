import openpyxl
import os
import tools
from mesh import TriangleMesh
from sequenceInfo import SequenceInfo


def write_sequences(sequence, E, D, M, rate_list, seq_ord, isFrame_I=False):
    for rate in rate_list:
        if isFrame_I == False:
            write_sequence(sequence, E[rate][0],
                           D[rate][0], M[rate][0], seq_ord)
        else:
            write_sequence_I(sequence, E[rate],
                             D[rate], M[rate], seq_ord)


def write_sequence_I(sequence, encoder_logs, decoder_logs, mmetric_logs, seq_ord):
    num_I_frames = len(encoder_logs)
    input_faces = []
    geometry_precision = []
    N0_frames = []
    frame_rate = []
    output_faces = []
    total_bitstream1 = []
    bmesh_intra = []
    bmesh_inter = []
    disp = []
    texture = []
    meta_data = []
    user_encoder_runtime = []
    user_decoder_runtime = []
    for i in range(num_I_frames):
        encoder_log = encoder_logs[i]
        decoder_log = decoder_logs[i]
        mmetric_log = mmetric_logs[i]
        s = seq_ord
        sxx, _ = tools.extract_variant(
            tools.extract_path(encoder_log, -2, -1), 's')
        rxx, r = tools.extract_variant(
            tools.extract_path(encoder_log, -2, -1), 'r')
        row = 5 + 5*s+r - 1
        row = str(row)

        input_faces.append(int(tools.extract_src_sth(
            mmetric_log, r'Triangles: (\d+)', isSum=True, onlyInput=True)))
        geometry_precision .append(int(tools.extract_src_sth(
            encoder_log, r'positionBitDepth\s+:\s+(\d+)')))
        N0_frames.append(int(tools.extract_src_sth(
            encoder_log, r'frameCount\s*:\s*(\d+)')))
        frame_rate.append(int(tools.extract_src_sth(
            encoder_log, r'framerate\s*:\s*(\d+)')))

        rate = '=IF(OR(ISBLANK(O'+row+'),ISBLANK(G'+row+')),"·",O' + \
            row+'*G'+row+'/F'+row+'/1000000)'
        total_bits_per_input_face = '=IF(ISBLANK(O' + \
            row+'),"·",O'+row+'/D'+row+')'
        total_bitstream = '=IF(OR(ISBLANK(M'+row+'),ISBLANK(S' + \
            row+'),ISBLANK(T'+row+')),".",M'+row+'+S'+row+'+T'+row+')'
        total_geo_bitstream = '=IF(OR(ISBLANK(Q'+row+'),ISBLANK(P' + \
            row+'),ISBLANK(R'+row+')),".",P'+row+'+Q'+row+'+R'+row+')'
        output_faces.append(int(tools.extract_src_sth(
            encoder_log, r'Sequence face count\s+(\d+)')))
        total_bitstream1.append(int(tools.extract_src_sth(
            encoder_log, r'Total:\s+\d+\s+B\s+(\d+)\s+b')))
        bmesh_intra.append(int(tools.extract_src_sth(
            encoder_log, r'Intra =\s+(\d+) B')) * 8)
        bmesh_inter.append(int(tools.extract_src_sth(
            encoder_log, r'Inter =\s+(\d+) B')) * 8)
        disp.append(int(tools.extract_src_sth(
            encoder_log, r'Displacement:\s+\d+\s+B\s+(\d+)\s+b')))
        texture .append(int(tools.extract_src_sth(
            encoder_log, r'Attribute:\s+\d+\s+B\s+(\d+)\s+b')))
        meta_data.append(int(tools.extract_src_sth(
            encoder_log, r'Metadata:\s+\d+\s+B\s+(\d+)\s+b')))
        psnr_D1 = float(tools.extract_src_sth(
            mmetric_log, r'Metric_results\s*:\s*([\d.]+)'))
        psnr_D2 = float(tools.extract_src_sth(
            mmetric_log, r'Metric_results\s*:\s+[\d.]+\s*([\d.]+)'))
        luma = float(tools.extract_src_sth(
            mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
        Chroma_Cb = float(tools.extract_src_sth(
            mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))

        Chroma_Cr = float(tools.extract_src_sth(
            mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
        image_base_geom = float(tools.extract_src_sth(
            mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
        image_base_luma = float(tools.extract_src_sth(
            mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
        user_encoder_runtime.append(int(float(tools.extract_src_sth(
            encoder_log, r'Sequence processing time\s+(\d+\.\d+) s'))))
        user_decoder_runtime.append(int(float(tools.extract_src_sth(
            decoder_log, r'Sequence processing time\s+(\d+\.\d+) s'))))

    sheet_name = 'C2 lossy RA I'
    sheet = workbook[sheet_name]

    sheet['A'+row] = sxx
    sheet['B'+row] = sequence
    sheet['c'+row] = rxx
    sheet['D'+row] = input_faces
    sheet['E'+row] = geometry_precision
    sheet['F'+row] = N0_frames
    sheet['G'+row] = frame_rate
    sheet['H'+row] = '=ROUND(J'+row+', 0)'

    sheet['J'+row] = rate
    sheet['K'+row] = total_bits_per_input_face
    sheet['L'+row] = total_bitstream
    sheet['M'+row] = total_geo_bitstream
    sheet['N'+row] = output_faces

    sheet['O'+row] = total_bitstream1
    sheet['P'+row] = bmesh_intra
    sheet['Q'+row] = bmesh_inter
    sheet['R'+row] = disp
    sheet['S'+row] = texture
    sheet['T'+row] = meta_data

    sheet['U'+row] = psnr_D1
    sheet['V'+row] = psnr_D2
    sheet['W'+row] = luma
    sheet['X'+row] = Chroma_Cb
    sheet['Y'+row] = Chroma_Cr

    sheet['Z'+row] = image_base_geom
    sheet['AA'+row] = image_base_luma

    sheet['AB'+row] = user_encoder_runtime
    sheet['AC'+row] = user_decoder_runtime


def write_sequence(sequence, encoder_log, decoder_log, mmetric_log, seq_ord):
    s = seq_ord
    sxx, _ = tools.extract_variant(
        tools.extract_path(encoder_log, -2, -1), 's')
    rxx, r = tools.extract_variant(
        tools.extract_path(encoder_log, -2, -1), 'r')
    row = 5 + 5*s+r - 1
    row = str(row)

    # src_mesh_path = tools.extract_src_sth(
    #     encoder_log, r'srcMesh\s*:\s*"(.*?)"')
    # frame_count, frame_start, frame_end = tools.count_files(src_mesh_path)

    # mesh0 = TriangleMesh()
    input_faces = int(tools.extract_src_sth(
        mmetric_log, r'Triangles: (\d+)', isSum=True, onlyInput=True))
    geometry_precision = int(tools.extract_src_sth(
        encoder_log, r'positionBitDepth\s+:\s+(\d+)'))
    N0_frames = int(tools.extract_src_sth(
        encoder_log, r'frameCount\s*:\s*(\d+)'))
    frame_rate = int(tools.extract_src_sth(
        encoder_log, r'framerate\s*:\s*(\d+)'))

    rate = '=IF(OR(ISBLANK(O'+row+'),ISBLANK(G'+row+')),"·",O' + \
        row+'*G'+row+'/F'+row+'/1000000)'
    total_bits_per_input_face = '=IF(ISBLANK(O'+row+'),"·",O'+row+'/D'+row+')'
    total_bitstream = '=IF(OR(ISBLANK(M'+row+'),ISBLANK(S' + \
        row+'),ISBLANK(T'+row+')),".",M'+row+'+S'+row+'+T'+row+')'
    total_geo_bitstream = '=IF(OR(ISBLANK(Q'+row+'),ISBLANK(P' + \
        row+'),ISBLANK(R'+row+')),".",P'+row+'+Q'+row+'+R'+row+')'
    output_faces = int(tools.extract_src_sth(
        encoder_log, r'Sequence face count\s+(\d+)'))
    total_bitstream1 = int(tools.extract_src_sth(
        encoder_log, r'Total:\s+\d+\s+B\s+(\d+)\s+b'))
    bmesh_intra = int(tools.extract_src_sth(
        encoder_log, r'Intra =\s+(\d+) B')) * 8
    bmesh_inter = int(tools.extract_src_sth(
        encoder_log, r'Inter =\s+(\d+) B')) * 8
    disp = int(tools.extract_src_sth(
        encoder_log, r'Displacement:\s+\d+\s+B\s+(\d+)\s+b'))
    texture = int(tools.extract_src_sth(
        encoder_log, r'Attribute:\s+\d+\s+B\s+(\d+)\s+b'))
    meta_data = int(tools.extract_src_sth(
        encoder_log, r'Metadata:\s+\d+\s+B\s+(\d+)\s+b'))
    psnr_D1 = float(tools.extract_src_sth(
        mmetric_log, r'Metric_results\s*:\s*([\d.]+)'))
    psnr_D2 = float(tools.extract_src_sth(
        mmetric_log, r'Metric_results\s*:\s+[\d.]+\s*([\d.]+)'))
    luma = float(tools.extract_src_sth(
        mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
    Chroma_Cb = float(tools.extract_src_sth(
        mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))

    Chroma_Cr = float(tools.extract_src_sth(
        mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
    image_base_geom = float(tools.extract_src_sth(
        mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
    image_base_luma = float(tools.extract_src_sth(
        mmetric_log, r'Metric_results\s*:\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*([\d.]+)'))
    user_encoder_runtime = int(float(tools.extract_src_sth(
        encoder_log, r'Sequence processing time\s+(\d+\.\d+) s')))
    user_decoder_runtime = int(float(tools.extract_src_sth(
        decoder_log, r'Sequence processing time\s+(\d+\.\d+) s')))

    sheet_name = 'C2 lossy RA'
    sheet = workbook[sheet_name]

    sheet['A'+row] = sxx
    sheet['B'+row] = sequence
    sheet['c'+row] = rxx
    sheet['D'+row] = input_faces
    sheet['E'+row] = geometry_precision
    sheet['F'+row] = N0_frames
    sheet['G'+row] = frame_rate
    sheet['H'+row] = '=ROUND(J'+row+', 0)'

    sheet['J'+row] = rate
    sheet['K'+row] = total_bits_per_input_face
    sheet['L'+row] = total_bitstream
    sheet['M'+row] = total_geo_bitstream
    sheet['N'+row] = output_faces

    sheet['O'+row] = total_bitstream1
    sheet['P'+row] = bmesh_intra
    sheet['Q'+row] = bmesh_inter
    sheet['R'+row] = disp
    sheet['S'+row] = texture
    sheet['T'+row] = meta_data

    sheet['U'+row] = psnr_D1
    sheet['V'+row] = psnr_D2
    sheet['W'+row] = luma
    sheet['X'+row] = Chroma_Cb
    sheet['Y'+row] = Chroma_Cr

    sheet['Z'+row] = image_base_geom
    sheet['AA'+row] = image_base_luma

    sheet['AB'+row] = user_encoder_runtime
    sheet['AC'+row] = user_decoder_runtime


dataSet_path = [
    "D:/mesh_data/voxelized/",
    "D:/mesh_data/v1/",
    "D:/mesh_data/v2/",
    "D:/mesh_data/v3/",
    "D:/mesh_data/v4/",]
dataList = [
    # 'basketball_player_voxelized/basketball_player_fr{:04d}_qp12_qt12.obj',
    # 'dancer_voxelized/dancer_fr{:04d}_qp12_qt12.obj',
    # 'football_voxelized/football_fr{:04d}_qp12_qt13.obj',
    'levi_voxelized/levi_fr{:04d}_qp12_qt13.obj',
    # 'longdress_voxelized/longdress_fr{:04d}_qp10_qt12.obj',
    'mitch_voxelized/mitch_fr{:04d}_qp12_qt13.obj',
    # 'soldier_voxelized/soldier_fr{:04d}_qp10_qt12.obj',
    'thomas_voxelized/thomas_fr{:04d}_qp12_qt13.obj',

    "cam1el-collapse/cam1el-collapse_fr{:04d}_qp12_qt13.obj",
    "cam2el-gallop/cam2el-gallop_fr{:04d}_qp12_qt13.obj",
    # "cam3el-poses/cam3el-poses_fr{:04d}_qp12_qt13.obj",
    # "catt-poses/catt-poses_fr{:04d}_qp12_qt13.obj",
    "ele2phant-gallop/ele2phant-gallop_fr{:04d}_qp12_qt13.obj",
    # "ele3phant-poses/ele3phant-poses_fr{:04d}_qp12_qt13.obj",
    # "flamingo-poses/flamingo-poses_fr{:04d}_qp12_qt13.obj",
    # "hor1se-collapse/hor1se-collapse_fr{:04d}_qp12_qt13.obj",#?????
    # "hor2se-gallop/hor2se-gallop_fr{:04d}_qp12_qt13.obj",#????
    # "hor3se-poses/hor3se-poses_fr{:04d}_qp12_qt13.obj",
    # "lion-poses/lion-poses_fr{:04d}_qp12_qt13.obj",

    "D_bouncing/D_bouncing_fr{:04d}_qp12_qt13.obj",
    "D_handstand/D_handstand_fr{:04d}_qp12_qt13.obj",
    "D_march/D_march_fr{:04d}_qp12_qt13.obj",
    "D_squat/D_squat_fr{:04d}_qp12_qt13.obj",
    "I_crane/I_crane_fr{:04d}_qp12_qt13.obj",
    "I_jumping/I_jumping_fr{:04d}_qp12_qt13.obj",
    "I_march/I_march_fr{:04d}_qp12_qt13.obj",
    "I_squat/I_squat_fr{:04d}_qp12_qt13.obj",
    "swing/swing_fr{:04d}_qp12_qt13.obj",
    "T_samba/T_samba_fr{:04d}_qp12_qt13.obj",

    # "chicken/chicken_fr{:04d}_qp12_qt13.obj",
    # "Cloth/Cloth_fr{:04d}_qp12_qt13.obj",  # no
    # "Dance/Dance_fr{:04d}_qp12_qt13.obj",  # no
    # "Jump/Jump_fr{:04d}_qp12_qt13.obj",
    # "mocapdance/mocapdance_fr{:04d}_qp12_qt13.obj",

    # "bunny/bunny_fr{:04d}_qp12_qt13.obj",
    # "james/james_fr{:04d}_qp12_qt13.obj",
    # "jessi/jessi_fr{:04d}_qp12_qt13.obj",
    # "nissan/nissan_fr{:04d}_qp12_qt13.obj",
    # "outt/outt_fr{:04d}_qp12_qt13.obj",
]

max_gof_size = 32

# 定义文件路径
file_path = os.path.dirname(os.path.abspath(
    __file__)) + "\\log_res\\F_all_v04.0.xlsm"


# 打开 Excel 文件
workbook = openpyxl.load_workbook(file_path)

# 选择特定的工作表
is_frames_I = [1]
rate_list = [1, 2, 3, 4, 5]

for is_frame_I in is_frames_I:
    seq_ord = 0
    for dataframe in dataList:
        input_path = "D:/mesh_data/voxelized/" + dataframe

        frame_count, frame_start, _ = tools.count_files(input_path)
        sequence = tools.extract_path(dataframe, -2, -1)
        print(str(seq_ord) + '\t:'+sequence)
        output_path = "log/" + sequence + \
            "_f" + str(frame_count) + "_gof.txt"

        seqInfo = SequenceInfo()
        seqInfo.load(frame_count, frame_start,
                     max_gof_size, output_path)
        if is_frame_I:
            E, D, M = seqInfo.getLogI(
                dataframe, frame_start, frame_count, rate_list)
            write_sequences(sequence, E, D, M, rate_list, seq_ord, is_frame_I)
            # seqInfo.encodeFrameI(dataframe, rate_list)
        else:
            E, D, M = seqInfo.getLogIP(
                dataframe, frame_start, frame_count, rate_list)
            write_sequences(sequence, E, D, M, rate_list, seq_ord, is_frame_I)
        seq_ord += 1

# encoder_log = "D:/mesh_code/mpeg-vmesh-tm/results/F032/s5c2r1_mitc/encoder.log"
# decoder_log = "D:/mesh_code/mpeg-vmesh-tm/results/F032/s5c2r1_mitc/decoder.log"
# mmetric_log = "D:/mesh_code/mpeg-vmesh-tm/results/F032/s5c2r1_mitc/metrics.log"


# sequence = "mitch_voxelized"


# 保存文件
workbook.save(file_path)

# 关闭工作簿
workbook.close()
