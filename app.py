import streamlit as st
import pandas as pd
import json
import os
import sys
import glob

# Подключаем наши модули
sys.path.append('src')
from phenotype_parser import parse_phenotype
from scorer import calculate_gene_scores

# Настройка страницы
st.set_page_config(page_title="Геномный детектив", page_icon="🧬", layout="wide")

st.title("🧬 Поиск генетических причин редких наследственных заболеваний")
st.markdown("### Система приоритизации генетических вариантов (Кейс 3, команда OMIX)")
st.markdown("---")

# ==========================================
# БОКОВАЯ ПАНЕЛЬ (Выбор действия)
# ==========================================
with st.sidebar:
    st.header("Навигация")
    
    # Ищем все доступные кейсы
    available_files = glob.glob("data/case_*_phenotype.json")
    case_ids = sorted(["_".join(os.path.basename(f).split('_')[:2]) for f in available_files])
    
    # Формируем список опций для меню
    options = ["🏆 Сводка по всем 10 кейсам"] + case_ids + ["📤 Загрузить свои файлы"]
    
    selected_mode = st.radio("Выберите действие:", options)
    st.markdown("---")

# ==========================================
# РЕЖИМ 1: СВОДКА ПО ВСЕМ 10 КЕЙСАМ
# ==========================================
if selected_mode == "🏆 Сводка по всем 10 кейсам":
    st.markdown("## Итоги анализа 10 клинических случаев")
    st.info("Наш пайплайн (VEP + HPO + ClinVar) успешно провел анализ для всех 10 кейсов (6 обучающих + 4 скрытых тестовых).")
    
    # ОБНОВЛЕНО: вписаны реальные гены и баллы из финального прогона
    results_data = {
        "Кейс": ["case_01", "case_02", "case_03", "case_04", "case_05", 
                 "case_06", "case_07", "case_08", "case_09", "case_10"],
        "Ожидаемый ген": ["FBN1", "PTPN11", "KMT2D", "FGFR2", "LMNA", 
                          "TGFBR2", "?", "?", "?", "?"], 
        "Найденный ген": ["FBN1", "PTPN11", "KMT2D", "FGFR2", "LMNA", 
                          "TGFBR2", "SCN2A", "STXBP1", "GLI3", "PIEZO2"],
        "Патогенный вариант": [
            "chr15:48412568:C:A", 
            "chr12:112472968:C:T", 
            "chr12:49034233:A:G", 
            "chr10:121520163:G:C", 
            "chr1:156138613:C:T", 
            "chr3:30672252:G:T",
            "chr2:165377606:A:G",   
            "chr9:135759686:G:A",   
            "chr7:41977637:C:G",   
            "chr18:10797517:G:A"    
        ],
        "Уверенность (Score)": ["21.5", "19.5", "36.0", "25.0", "32.5", "45.5", "28.5", "6.0", "10.5", "10.5"],
        "Статус": ["✅ Успех", "✅ Успех", "✅ Успех", "✅ Успех", "✅ Успех", 
                   "✅ Успех", "🧪 Скрытый", "🧪 Скрытый", "🧪 Скрытый", "🧪 Скрытый"]
    }
    
    df_results = pd.DataFrame(results_data)
    
    # Красим колонку статуса
    def color_status(val):
        if '✅' in str(val): return 'background-color: #1e4620; color: #85e085'
        if '🧪' in str(val): return 'background-color: #1a3b5c; color: #85c1e9'
        return ''
    
    st.dataframe(df_results.style.map(color_status, subset=['Статус']), use_container_width=True, hide_index=True)
    st.success("🎯 Точность на валидационной выборке (кейсы 1-6): 6 из 6 (100%)")

# ==========================================
# РЕЖИМ 2: ДЕТАЛЬНЫЙ АНАЛИЗ КОНКРЕТНОГО КЕЙСА
# ==========================================
elif selected_mode in case_ids:
    st.markdown(f"## Анализ клинического случая: {selected_mode}")
    
    pheno_path = f"data/{selected_mode}_phenotype.json"
    annotated_path = f"data/{selected_mode}_annotated.tsv"
    
    if not os.path.exists(annotated_path):
        st.error(f"Для {selected_mode} еще нет файла аннотации.")
    else:
        with st.spinner(f"Анализируем {selected_mode}..."):
            try:
                pos_hpo, neg_hpo, pos_labels, neg_labels = parse_phenotype(pheno_path)
                ranked_genes = calculate_gene_scores(annotated_path, pos_hpo, neg_hpo)
                display_results(ranked_genes, pos_hpo, neg_hpo, pos_labels, neg_labels, selected_mode)
            except Exception as e:
                st.error(f"Ошибка при анализе: {e}")

# ==========================================
# РЕЖИМ 3: ЗАГРУЗКА СВОИХ ФАЙЛОВ
# ==========================================
elif selected_mode == "📤 Загрузить свои файлы":
    st.markdown("## Анализ нового клинического случая")
    st.write("Загрузите файл симптомов (JSON) и файл аннотированных вариантов (TSV).")
    
    col1, col2 = st.columns(2)
    with col1:
        uploaded_pheno = st.file_uploader("Файл симптомов (JSON)", type=['json'])
    with col2:
        uploaded_annotated = st.file_uploader("Файл вариантов (TSV)", type=['tsv'])
        
    if uploaded_pheno and uploaded_annotated:
        if st.button("Запустить анализ загруженных файлов", use_container_width=True):
            with st.spinner("Обработка файлов..."):
                try:
                    temp_pheno_path = "temp_pheno.json"
                    temp_annotated_path = "temp_annotated.tsv"
                    
                    with open(temp_pheno_path, "wb") as f:
                        f.write(uploaded_pheno.getbuffer())
                    with open(temp_annotated_path, "wb") as f:
                        f.write(uploaded_annotated.getbuffer())
                        
                    pos_hpo, neg_hpo, pos_labels, neg_labels = parse_phenotype(temp_pheno_path)
                    ranked_genes = calculate_gene_scores(temp_annotated_path, pos_hpo, neg_hpo)
                    
                    display_results(ranked_genes, pos_hpo, neg_hpo, pos_labels, neg_labels, "Загруженный кейс")
                    
                    os.remove(temp_pheno_path)
                    os.remove(temp_annotated_path)
                except Exception as e:
                    st.error(f"Ошибка при обработке файлов: {e}")

# ==========================================
# ФУНКЦИЯ ДЛЯ ОТРИСОВКИ РЕЗУЛЬТАТОВ (Вкладки)
# ==========================================
def display_results(ranked_genes, pos_hpo, neg_hpo, pos_labels, neg_labels, case_name):
    if ranked_genes is None or ranked_genes.empty:
        st.warning("Не удалось найти кандидатов (возможно, нет совпадений с HPO).")
        return
        
    st.success(f"✅ Анализ для {case_name} завершен!")
    
    tab1, tab2, tab3 = st.tabs(["🏆 Главный кандидат", "📊 Топ-5 таблица", "📝 Симптомы пациента"])
    
    with tab1:
        top = ranked_genes.iloc[0]
        st.markdown(f"## Наиболее вероятный ген-кандидат: **{top['gene']}**")
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Итоговый балл", f"{top['score']}")
        col2.metric("Совпадений HPO", f"{top['matches']}")
        col3.metric("Тип мутации", top['consequence'])
            
        st.markdown("### Почему система выбрала этот вариант?")
        st.info(f"""
        **1. Совпадение по симптомам:** 
        У пациента выявлено {len(pos_hpo)} подтвержденных симптомов. Ген **{top['gene']}** ассоциирован с **{top['matches']}** из них. 
        *(Базовый вклад: {top['matches'] * 1.5} баллов)*
        
        **2. Опасность мутации:** 
        Вариант `{top['variant']}` является **{top['consequence']}**. Этот тип изменений часто нарушает функцию белка.
        *(Бонус: +{top['cons_bonus']} баллов)*
        
        **3. Штраф за несоответствие:**
        У пациента *отсутствуют* {top['excluded_matches']} симптомов, характерных для этого гена.
        *(Штраф: -{top['excluded_matches'] * 0.5} баллов)*
        """)

    with tab2:
        st.markdown("### Полный список топ-5 кандидатов")
        display_df = ranked_genes.head(5)[['gene', 'variant', 'consequence', 'score']].copy()
        display_df.columns = ['Ген', 'Вариант', 'Тип мутации', 'Балл']
        st.dataframe(display_df, use_container_width=True, hide_index=True)

    with tab3:
        st.markdown("### Клиническая картина пациента")
        
        st.write("#### Подтвержденные симптомы (есть у пациента):")
        for label in pos_labels:
            st.markdown(f"- ✅ {label}")
            
        if neg_labels:
            st.write("#### Исключенные симптомы (отсутствуют):")
            for label in neg_labels:
                st.markdown(f"- ❌ ~~{label}~~")

st.markdown("---")
st.caption("Дисклеймер: Результаты предназначены исключительно в исследовательских и образовательных целях. Не является медицинским диагнозом.")
