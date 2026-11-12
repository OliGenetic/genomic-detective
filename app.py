import streamlit as st
import pandas as pd
import json
import os
import sys
import glob

sys.path.append('src')
from phenotype_parser import parse_phenotype
from scorer import calculate_gene_scores

st.set_page_config(page_title="Поиск генетических причин редких наследственных заболеваний", page_icon="🧬", layout="wide")

st.title("🧬 Поиск генетических причин редких наследственных заболеваний")
st.markdown("### Система приоритизации генетических вариантов для поиска редких заболеваний")
st.markdown("---")

with st.sidebar:
    st.header("Выбор кейса для анализа")
    available_files = glob.glob("data/case_*_phenotype.json")
    case_ids = sorted(["_".join(os.path.basename(f).split('_')[:2]) for f in available_files])
    
    selected_case = st.selectbox("Выберите пациента (кейс):", case_ids)
    run_button = st.button(" Запустить анализ", use_container_width=True)

if run_button:
    if not selected_case:
        st.error(" Пожалуйста, выберите кейс.")
    else:
        pheno_path = f"data/{selected_case}_phenotype.json"
        annotated_path = f"data/{selected_case}_annotated.tsv"
        
        if not os.path.exists(annotated_path):
            st.error(f"Для {selected_case} еще нет аннотации.")
        else:
            with st.spinner(f"Анализируем {selected_case}..."):
                try:
                    # Получаем теперь 4 переменные!
                    pos_hpo, neg_hpo, pos_labels, neg_labels = parse_phenotype(pheno_path)
                    ranked_genes = calculate_gene_scores(annotated_path, pos_hpo, neg_hpo)
                    
                    st.success(f" Анализ для {selected_case} завершен!")
                    
                    tab1, tab2, tab3 = st.tabs(["🏆 Главный кандидат", " Топ-5 таблица", "📝 Симптомы пациента"])
                    
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
                        # Выводим списком с галочками
                        for label in pos_labels:
                            st.markdown(f"- {label}")
                            
                        if neg_labels:
                            st.write("#### Исключенные симптомы (отсутствуют):")
                            for label in neg_labels:
                                st.markdown(f"- ~~{label}~~")

                except Exception as e:
                    st.error(f"Ошибка: {e}")

st.markdown("---")
st.caption("Дисклеймер: Результаты предназначены исключительно в исследовательских и образовательных целях.")