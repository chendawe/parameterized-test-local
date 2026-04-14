# pip install chromadb langchain langchain-community sentence_transformers

import os
from langchain.text_splitter import CharacterTextSplitter
from langchain_community.document_loaders import TextLoader
from langchain_community.vectorstores import Chroma
# 假设你有一个存放文本文件的目录，比如 ./books
books_dir = "./books"
persist_directory = "./my_chroma_db"  # Chroma数据持久化目录

# 1. 加载文档
documents = []
for file in os.listdir(books_dir):
    if file.endswith(".txt"):
        loader = TextLoader(os.path.join(books_dir, file))
        documents.extend(loader.load())

# 2. 文本分块 (Chunking)
text_splitter = CharacterTextSplitter(chunk_size=500, chunk_overlap=50)
docs = text_splitter.split_documents(documents)

# 3. 创建向量库并持久化
#    这里我们没有指定embedding函数，LangChain会默认使用sentence_transformers的一个模型
vectordb = Chroma.from_documents(
    documents=docs, 
    embedding=None,  # 留空则使用默认
    persist_directory=persist_directory
)
vectordb.persist()  # 确保数据写入磁盘
print(f"向量库创建完成，共 {len(docs)} 个文本块，已存储在 {persist_directory}")






from langchain_community.vectorstores import Chroma
from langchain_openai import ChatOpenAI  # 或其他LLM
from langchain.chains import RetrievalQA

# 1. 加载已存在的向量库（需要指定相同的embedding函数）
#    这里我们假设使用的是默认embedding，加载时必须保持一致
vectordb = Chroma(persist_directory=persist_directory, embedding_function=None)

# 2. 将向量库转换为检索器
retriever = vectordb.as_retriever(search_kwargs={"k": 3})  # 检索最相关的3个片段

# 3. 创建问答链
llm = ChatOpenAI(model="gpt-3.5-turbo", temperature=0)  # 假设你有OpenAI key
qa_chain = RetrievalQA.from_chain_type(
    llm=llm, 
    chain_type="stuff",  # 简单的将检索到的文档“塞”进提示词
    retriever=retriever
)

# 4. 提问并获取答案
query = "这本书里提到了哪些主要观点？"
response = qa_chain.invoke(query)
print("回答：", response['result'])

# 你还可以直接查看检索到的原始文档内容
retrieved_docs = retriever.invoke(query)
print("\n--- 检索到的相关文档片段 ---")
for i, doc in enumerate(retrieved_docs):
    print(f"片段 {i+1}:\n{doc.page_content}\n")